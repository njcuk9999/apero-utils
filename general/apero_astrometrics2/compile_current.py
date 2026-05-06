#!/usr/bin/env python3
"""Seed astrometrics YAML files from a CSV catalog."""

from __future__ import annotations

import os
import argparse
import csv
import re
from pathlib import Path
from typing import Any

import yaml


DEFAULT_INPUT_CANDIDATES = ("main.csv", "main_list.csv")
DEFAULT_OUTPUT_DIR = "astrometrics"


def count_data_rows(csv_path: Path) -> int:
    """Count CSV data rows (excluding header) for progress reporting."""
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        return sum(1 for _ in handle) - 1


def normalize_value(value: Any) -> Any:
    """Convert empty-like values to None so YAML writes null."""
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"none", "null", "nan"}:
        return None
    return text


def to_float_or_none(value: Any) -> float | None:
    value = normalize_value(value)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_list_field(value: Any) -> list[str] | None:
    value = normalize_value(value)
    if value is None:
        return None
    items = [item.strip() for item in str(value).split("|") if item.strip()]
    return items or None


def classify_apero_class(keywords: list[str] | None) -> str:
    """Default STAR, overridden by FIELD / SOLAR_SYSTEM in KEYWORDS."""
    if not keywords:
        return "STAR"
    upper = {item.upper() for item in keywords}
    if "SOLAR_SYSTEM" in upper:
        return "SOLAR_SYSTEM"
    if "FIELD" in upper:
        return "FIELD"
    return "STAR"


def sanitize_filename(objname: str) -> str:
    """Keep filenames safe across filesystems."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", objname).strip("._") or "unnamed"


def scalar_with_source(value: Any, source: Any, units: str | None = None) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "value": to_float_or_none(value),
        "source": normalize_value(source),
    }
    if units is not None:
        payload["units"] = units
    return payload


def text_with_source(value: Any, source: Any) -> dict[str, Any]:
    return {
        "value": normalize_value(value),
        "source": normalize_value(source),
    }


def build_record(row: dict[str, str]) -> dict[str, Any]:
    keywords = parse_list_field(row.get("KEYWORDS"))

    return {
        "APERO_NAME": normalize_value(row.get("OBJNAME")),
        "ORIGINAL_NAME": normalize_value(row.get("ORIGINAL_NAME")),
        "SIMBAD_NAME": None,
        "APERO_CLASS": classify_apero_class(keywords),
        "RA": scalar_with_source(row.get("RA_DEG"), row.get("RA_SOURCE"), "deg"),
        "DEC": scalar_with_source(row.get("DEC_DEG"), row.get("DEC_SOURCE"), "deg"),
        "EPOCH": to_float_or_none(row.get("EPOCH")),
        "PMRA": scalar_with_source(row.get("PMRA"), row.get("PMRA_SOURCE"), "mas/yr"),
        "PMDE": scalar_with_source(row.get("PMDE"), row.get("PMDE_SOURCE"), "mas/yr"),
        "PLX": scalar_with_source(row.get("PLX"), row.get("PLX_SOURCE"), "mas"),
        "RV": scalar_with_source(row.get("RV"), row.get("RV_SOURCE"), "km/s"),
        "TEFF": scalar_with_source(row.get("TEFF"), row.get("TEFF_SOURCE"), "K"),
        "SPT": text_with_source(row.get("SP_TYPE"), row.get("SP_SOURCE")),
        "VSINI": {
            "value": to_float_or_none(row.get("VSINI")),
            "err": to_float_or_none(row.get("VSINI_ERR")),
            "source": normalize_value(row.get("VSINI_SOURCE")),
            "units": "km/s",
        },
        "G_MAG": text_with_source(row.get("GMAG"), row.get("GMAG_SOURCE")),
        "GBP_MAG": text_with_source(None, None),
        "GRP_MAG": text_with_source(None, None),
        "J_MAG": text_with_source(row.get("JMAG"), row.get("JMAG_SOURCE")),
        "H_MAG": text_with_source(row.get("HMAG"), row.get("HMAG_SOURCE")),
        "KS_MAG": text_with_source(row.get("KMAG"), row.get("KMAG_SOURCE")),
        "W1_MAG": text_with_source(None, None),
        "W2_MAG": text_with_source(None, None),
        "W3_MAG": text_with_source(None, None),
        "KEYWORDS": keywords,
        "ALIASES": parse_list_field(row.get("ALIASES")),
        "NOTES": normalize_value(row.get("NOTES")),
    }


def find_input_csv(script_dir: Path, explicit_input: str | None) -> Path:
    if explicit_input:
        return Path(explicit_input).expanduser().resolve()
    for candidate in DEFAULT_INPUT_CANDIDATES:
        path = script_dir / candidate
        if path.exists():
            return path
    return (script_dir / DEFAULT_INPUT_CANDIDATES[0]).resolve()


def write_yaml(path: Path, payload: dict[str, Any]) -> None:
    text = yaml.safe_dump(payload, sort_keys=False, allow_unicode=False)
    path.write_text(text, encoding="utf-8")


def run(
    input_csv: Path,
    output_dir: Path,
    overwrite: bool,
    limit: int | None,
    progress_every: int,
) -> tuple[int, int]:
    output_dir.mkdir(parents=True, exist_ok=True)

    written = 0
    skipped = 0

    total_rows = count_data_rows(input_csv)
    if limit is not None:
        total_rows = min(total_rows, limit)

    with input_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for index, row in enumerate(reader, start=1):
            if limit is not None and index > limit:
                break

            objname = normalize_value(row.get("OBJNAME"))
            if objname is None:
                skipped += 1
                continue

            outfile = output_dir / f"{sanitize_filename(str(objname))}.yaml"
            if outfile.exists() and not overwrite:
                skipped += 1
                continue

            payload = build_record(row)
            write_yaml(outfile, payload)
            written += 1

            if progress_every > 0 and (index % progress_every == 0 or index == total_rows):
                print(
                    f"Processed {index}/{total_rows} rows "
                    f"(written={written}, skipped={skipped})",
                    end="\r",
                    flush=True,
                )

    if progress_every > 0:
        print()

    return written, skipped


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default=None, help="Path to CSV (default: auto main.csv/main_list.csv)")
    parser.add_argument(
        "--output-dir",
        default=str(script_dir / DEFAULT_OUTPUT_DIR),
        help="Output directory for OBJNAME.yaml files",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing YAML files")
    parser.add_argument("--limit", type=int, default=None, help="Only process first N rows (debug)")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=500,
        help="Print progress every N rows (set 0 to disable)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent

    input_csv = find_input_csv(script_dir, args.input)
    output_dir = Path(args.output_dir).expanduser().resolve()

    written, skipped = run(
        input_csv=input_csv,
        output_dir=output_dir,
        overwrite=args.overwrite,
        limit=args.limit,
        progress_every=args.progress_every,
    )
    print(f"Input CSV: {input_csv}")
    print(f"Output dir: {output_dir}")
    print(f"Wrote {written} YAML files (skipped {skipped}).")


if __name__ == "__main__":
    try:
        _ = __file__
    except NameError:
        __file__ = os.path.join(os.getcwd(), 'compile_current.py')
    main()

