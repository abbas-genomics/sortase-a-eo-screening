"""
extract.py - PubChem extraction pipeline
=======================================
Copyright (c) Abbas Aliyu. All rights reserved.
Developed by Abbas Aliyu Abbas
"""

import argparse
import csv
import logging
import re
import shutil
import time
import unicodedata
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote

import requests


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CACTUS_BASE = "https://cactus.nci.nih.gov/chemical/structure"

# Properties to request from PubChem PUG REST API.
# Note: PubChem returns "ConnectivitySMILES" for CanonicalSMILES requests
# and "SMILES" for IsomericSMILES requests — we remap after fetching.
PROPERTIES_REQUEST = [
    "MolecularFormula",
    "MolecularWeight",
    "IUPACName",
    "XLogP",
    "TPSA",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
    "IsomericSMILES",
    "CanonicalSMILES",
    "InChIKey",
]

# Map from PubChem JSON response key -> column name used in output CSVs.
PROPERTY_KEY_MAP = {
    "MolecularFormula": "MolecularFormula",
    "MolecularWeight": "MolecularWeight",
    "IUPACName": "IUPACName",
    "XLogP": "XLogP",
    "TPSA": "TPSA",
    "HBondDonorCount": "HBondDonorCount",
    "HBondAcceptorCount": "HBondAcceptorCount",
    "RotatableBondCount": "RotatableBondCount",
    "SMILES": "IsomericSMILES",          # PubChem returns SMILES for IsomericSMILES
    "ConnectivitySMILES": "CanonicalSMILES",  # PubChem returns ConnectivitySMILES for CanonicalSMILES
    "InChIKey": "InChIKey",
}

# Ordered output column names for properties CSV.
PROPERTIES_COLUMNS = [
    "MolecularFormula",
    "MolecularWeight",
    "IUPACName",
    "XLogP",
    "TPSA",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
    "IsomericSMILES",
    "CanonicalSMILES",
    "InChIKey",
]


def setup_logger(log_file: Path) -> logging.Logger:
    logger = logging.getLogger("pubchem_pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger


def sanitize_filename(name: str) -> str:
    replacements = {
        "ß": "beta",
        "β": "beta",
        "α": "alpha",
        "δ": "delta",
        "γ": "gamma",
    }
    for src, dst in replacements.items():
        name = name.replace(src, dst)

    cleaned = unicodedata.normalize("NFKD", name)
    cleaned = cleaned.encode("ascii", "ignore").decode("ascii")
    cleaned = cleaned.strip().replace(" ", "_")
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", cleaned)
    return cleaned.strip("._") or "compound"


def unique_sdf_path(base_dir: Path, base_name: str) -> Path:
    path = base_dir / f"{base_name}.sdf"
    if not path.exists():
        return path

    counter = 2
    while True:
        candidate = base_dir / f"{base_name}_{counter}.sdf"
        if not candidate.exists():
            return candidate
        counter += 1


def normalize_name_variants(name: str) -> List[str]:
    stripped = name.strip()
    variants = [stripped]
    replacements = {
        "ß": "beta",
        "β": "beta",
        "α": "alpha",
        "δ": "delta",
        "γ": "gamma",
    }

    replaced = stripped
    for src, dst in replacements.items():
        replaced = replaced.replace(src, dst)
    variants.append(replaced)

    # Handle names that include trailing acronyms in parentheses, e.g. "...(DEHP)".
    trailing_acronym_removed = re.sub(r"\(([A-Za-z0-9-]{2,10})\)$", "", replaced).strip()
    if trailing_acronym_removed and trailing_acronym_removed != replaced:
        variants.append(trailing_acronym_removed)

    # Also try with spaces around parentheses removed and punctuation softened.
    variants.append(replaced.replace("(", " ").replace(")", " "))
    variants.append(replaced.replace(",", " "))

    # Frequent typo/variant: methylen vs methylene.
    variants.append(replaced.replace("methylen", "methylene"))

    if "-" in replaced:
        variants.append(replaced.replace("-", " "))

    # PubChem occasionally stores synonyms with a trailing marker '#'.
    hashtag_variants = [f"{v} #" for v in variants if v]
    variants.extend(hashtag_variants)

    # Normalize whitespace and deduplicate while preserving order.
    cleaned_variants = []
    for variant in variants:
        v = re.sub(r"\s+", " ", variant).strip()
        if v:
            cleaned_variants.append(v)
    return list(dict.fromkeys(cleaned_variants))


def request_json(session: requests.Session, url: str, logger: logging.Logger) -> Optional[Dict[str, Any]]:
    try:
        response = session.get(url, timeout=30)
        if response.status_code != 200:
            return None
        return response.json()
    except requests.RequestException as exc:
        logger.warning("Request failed: %s | %s", url, exc)
        return None


def cid_from_name(session: requests.Session, candidate_name: str, logger: logging.Logger) -> Optional[int]:
    encoded = quote(candidate_name)
    url = f"{PUBCHEM_BASE}/compound/name/{encoded}/cids/JSON"
    data = request_json(session, url, logger)
    if not data:
        return None
    cids = data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return None
    return int(cids[0])


def search_pubchem_cid(session: requests.Session, compound_name: str, logger: logging.Logger) -> Tuple[Optional[int], Optional[str]]:
    for candidate in normalize_name_variants(compound_name):
        cid = cid_from_name(session, candidate, logger)
        if cid:
            return cid, candidate

    return None, None


def fetch_pubchem_properties(session: requests.Session, cid: int, logger: logging.Logger) -> Dict[str, Any]:
    prop_list = ",".join(PROPERTIES_REQUEST)
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/{prop_list}/JSON"
    data = request_json(session, url, logger)
    if not data:
        return {}

    rows = data.get("PropertyTable", {}).get("Properties", [])
    raw = rows[0] if rows else {}
    # Remap response keys to canonical output column names.
    return {PROPERTY_KEY_MAP.get(k, k): v for k, v in raw.items()}


def download_pubchem_sdf(
    session: requests.Session,
    cid: int,
    output_file: Path,
    logger: logging.Logger,
) -> Tuple[bool, Optional[str], Optional[str]]:
    for record_type in ("3d", "2d"):
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type={record_type}"
        try:
            response = session.get(url, timeout=60)
            if response.status_code == 200 and response.content:
                output_file.write_bytes(response.content)
                return True, record_type.upper(), "PubChem"
            if response.status_code not in (404, 400):
                logger.warning("SDF download returned status %s for CID %s", response.status_code, cid)
        except requests.RequestException as exc:
            logger.warning("PubChem SDF request failed for CID %s: %s", cid, exc)
    return False, None, None


def download_cactus_sdf(
    session: requests.Session,
    compound_name: str,
    output_file: Path,
    logger: logging.Logger,
) -> Tuple[bool, Optional[str], Optional[str]]:
    encoded = quote(compound_name)
    url = f"{CACTUS_BASE}/{encoded}/sdf"
    try:
        response = session.get(url, timeout=60)
        if response.status_code == 200 and response.text and "M  END" in response.text:
            output_file.write_text(response.text, encoding="utf-8")
            return True, "2D", "NCI_CACTUS"
    except requests.RequestException as exc:
        logger.warning("CACTUS SDF request failed for %s: %s", compound_name, exc)
    return False, None, None


def read_compounds(input_path: Path) -> List[str]:
    suffix = input_path.suffix.lower()

    if suffix in (".txt", ""):
        names = [line.strip() for line in input_path.read_text(encoding="utf-8").splitlines()]
        return [n for n in names if n]

    if suffix == ".csv":
        names: List[str] = []
        with input_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            for row in reader:
                if not row:
                    continue
                candidate = row[0].strip()
                if candidate:
                    names.append(candidate)
        return names

    if suffix in (".xlsx", ".xls"):
        try:
            import pandas as pd
        except ImportError as exc:
            raise RuntimeError("Reading Excel input requires pandas and openpyxl.") from exc

        df = pd.read_excel(input_path)
        first_col = df.columns[0]
        return [str(v).strip() for v in df[first_col].dropna().tolist() if str(v).strip()]

    raise ValueError(f"Unsupported input extension: {suffix}")


def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_pipeline(input_file: Path, output_dir: Path, delay: float = 0.25) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    sdf_dir = output_dir / "sdf_files"
    numbered_dir = output_dir / "numbered_sdf"
    sdf_dir.mkdir(parents=True, exist_ok=True)
    numbered_dir.mkdir(parents=True, exist_ok=True)

    log_file = output_dir / "pipeline.log"
    status_csv = output_dir / "download_status.csv"
    props_csv = output_dir / "compound_properties.csv"
    rename_log_csv = numbered_dir / "renaming_log.csv"
    smiles_file = numbered_dir / "smiles.csv"
    logger = setup_logger(log_file)

    compounds = read_compounds(input_file)
    if not compounds:
        logger.error("No compounds found in input file: %s", input_file)
        return 1

    logger.info("Starting pipeline for %d compounds", len(compounds))
    session = requests.Session()
    session.headers.update({"User-Agent": "pubchem-extract-pipeline/1.0"})

    status_rows: List[Dict[str, Any]] = []
    prop_rows: List[Dict[str, Any]] = []
    rename_rows: List[Dict[str, Any]] = []
    smiles_rows: List[Dict[str, Any]] = []

    success = 0
    failed = 0
    numbered_counter = 1

    for idx, raw_name in enumerate(compounds, start=1):
        compound = raw_name.strip()
        safe_name = sanitize_filename(compound)
        logger.info("[%d/%d] Processing: %s", idx, len(compounds), compound)

        cid, matched_name = search_pubchem_cid(session, compound, logger)
        if cid:
            logger.info("Resolved CID %s using name variant: %s", cid, matched_name)
        else:
            logger.warning("No PubChem CID found for: %s", compound)

        properties: Dict[str, Any] = {}
        if cid:
            properties = fetch_pubchem_properties(session, cid, logger)

        cid_prefix = str(cid) if cid else "NA"
        base_filename = f"{cid_prefix}_{safe_name}"
        sdf_path = unique_sdf_path(sdf_dir, base_filename)
        downloaded, structure_type, source = (False, None, None)

        if cid:
            downloaded, structure_type, source = download_pubchem_sdf(session, cid, sdf_path, logger)

        if not downloaded:
            downloaded, structure_type, source = download_cactus_sdf(session, compound, sdf_path, logger)

        status = "DONE" if downloaded else "ERROR"
        if downloaded:
            success += 1
            logger.info(
                "Downloaded %s structure for %s from %s -> %s",
                structure_type,
                compound,
                source,
                sdf_path.name,
            )

            numbered_name = f"{numbered_counter}.sdf"
            numbered_path = numbered_dir / numbered_name
            shutil.copy2(sdf_path, numbered_path)
            canonical_smiles = properties.get("CanonicalSMILES", "")
            isomeric_smiles = properties.get("IsomericSMILES", "")
            rename_rows.append(
                {
                    "index": numbered_counter,
                    "numbered_file": numbered_name,
                    "original_file": sdf_path.name,
                    "cid": cid or "",
                    "input_name": compound,
                    "structure_type": structure_type or "",
                    "source": source or "",
                    "CanonicalSMILES": canonical_smiles,
                    "IsomericSMILES": isomeric_smiles,
                }
            )
            smiles_rows.append(
                {
                    "index": numbered_counter,
                    "numbered_file": numbered_name,
                    "cid": cid or "",
                    "input_name": compound,
                    "CanonicalSMILES": canonical_smiles,
                    "IsomericSMILES": isomeric_smiles,
                    "structure_type": structure_type or "",
                }
            )
            numbered_counter += 1
        else:
            failed += 1
            logger.error("Could not download SDF for: %s", compound)

        status_rows.append(
            {
                "input_name": compound,
                "matched_name": matched_name or "",
                "cid": cid or "",
                "status": status,
                "structure_type": structure_type or "",
                "source": source or "",
                "sdf_file": str(sdf_path.name) if downloaded else "",
            }
        )

        prop_row: Dict[str, Any] = {
            "input_name": compound,
            "matched_name": matched_name or "",
            "cid": cid or "",
        }
        for key in PROPERTIES_COLUMNS:
            prop_row[key] = properties.get(key, "")
        prop_rows.append(prop_row)

        time.sleep(delay)

    write_csv(
        status_csv,
        status_rows,
        ["input_name", "matched_name", "cid", "status", "structure_type", "source", "sdf_file"],
    )
    write_csv(props_csv, prop_rows, ["input_name", "matched_name", "cid", *PROPERTIES_COLUMNS])
    write_csv(
        rename_log_csv,
        rename_rows,
        ["index", "numbered_file", "original_file", "cid", "input_name", "structure_type", "source", "CanonicalSMILES", "IsomericSMILES"],
    )
    write_csv(
        smiles_file,
        smiles_rows,
        ["index", "numbered_file", "cid", "input_name", "CanonicalSMILES", "IsomericSMILES", "structure_type"],
    )

    logger.info("Pipeline completed. DONE=%d | ERROR=%d | TOTAL=%d", success, failed, len(compounds))

    print("\n=== PIPELINE SUMMARY ===")
    print(f"Input file      : {input_file}")
    print(f"Output dir      : {output_dir}")
    print(f"Done compounds  : {success}")
    print(f"Error compounds : {failed}")
    print(f"Status CSV      : {status_csv}")
    print(f"Properties CSV  : {props_csv}")
    print(f"Numbered folder : {numbered_dir}")
    print(f"Renaming log    : {rename_log_csv}")
    print(f"SMILES file     : {smiles_file}")
    print(f"Log file        : {log_file}")

    return 0 if failed == 0 else 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Search compounds by name, download SDF (3D preferred, 2D fallback), "
            "and export status + properties CSV files."
        )
    )
    parser.add_argument(
        "--input",
        default="compounds_name",
        help="Input file with compound names (txt, csv, xlsx). Default: compounds_name",
    )
    parser.add_argument(
        "--output-dir",
        default="output",
        help="Directory for SDF files, logs, and CSV outputs. Default: output",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.25,
        help="Delay between requests in seconds. Default: 0.25",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output_dir)

    if not input_path.exists():
        print(f"Input file not found: {input_path}")
        raise SystemExit(1)

    raise SystemExit(run_pipeline(input_path, output_path, delay=args.delay))