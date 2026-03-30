"""Startup data bootstrap utilities.

Downloads required input artifacts on first run and skips work when data already exists.
"""

from __future__ import annotations

import logging
import os
import shutil
import tarfile
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict
from urllib.error import URLError
from urllib.request import urlopen

LOGGER = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent
SOURCES_FILE = PROJECT_ROOT / "data_sources.txt"
RESULTS_FILENAME = "cis_trans_results_table.csv"
GENE_DATA_ARCHIVE = "gene_count_data.tar.gz"
GENE_DATA_DIRNAME = "gene_count_data"


@dataclass
class BootstrapResult:
    downloaded: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    details: list[str] = field(default_factory=list)
    error: str | None = None


def get_data_root() -> Path:
    """Resolve data root from DATA_ROOT env var; defaults to project root."""
    configured = os.getenv("DATA_ROOT", ".").strip() or "."
    path = Path(configured)
    if not path.is_absolute():
        path = (PROJECT_ROOT / path).resolve()
    return path


def _load_sources(path: Path = SOURCES_FILE) -> Dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing data source file: {path}")

    sources: Dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        if "|" not in line:
            raise ValueError(f"Invalid data source entry: {raw_line}")

        artifact, url = [part.strip() for part in line.split("|", 1)]
        if not artifact or not url:
            raise ValueError(f"Invalid data source entry: {raw_line}")
        sources[artifact] = url

    required = {RESULTS_FILENAME, GENE_DATA_ARCHIVE}
    missing = required.difference(sources)
    if missing:
        missing_list = ", ".join(sorted(missing))
        raise ValueError(f"data_sources.txt is missing required entries: {missing_list}")
    return sources


def _download_to_path(url: str, destination: Path, timeout: int = 120) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=timeout) as response, destination.open("wb") as out_file:
        shutil.copyfileobj(response, out_file)


def _is_safe_member(member_name: str, extract_root: Path) -> bool:
    destination = (extract_root / member_name).resolve()
    root = extract_root.resolve()
    return destination == root or root in destination.parents


def _safe_extract_tar(archive_path: Path, extract_root: Path) -> None:
    with tarfile.open(archive_path, mode="r:gz") as tar:
        members = tar.getmembers()
        for member in members:
            if not _is_safe_member(member.name, extract_root):
                raise ValueError(f"Unsafe path in archive member: {member.name}")
        tar.extractall(path=extract_root, members=members)


def _gene_data_ready(data_root: Path) -> bool:
    subtype_dir = data_root / GENE_DATA_DIRNAME / "subtype"
    if not subtype_dir.is_dir():
        return False

    for tissue_dir in subtype_dir.iterdir():
        if not tissue_dir.is_dir():
            continue
        has_counts = False
        has_meta = False
        for child in tissue_dir.iterdir():
            name = child.name
            if name.endswith("_xgener_input_dataframe_FILTERED.csv"):
                has_counts = True
            elif name.endswith("_xgener_input_metadata_FILTERED.csv"):
                has_meta = True
            if has_counts and has_meta:
                return True
    return False


def _find_extracted_gene_data_dir(extract_root: Path) -> Path:
    direct = extract_root / GENE_DATA_DIRNAME
    if direct.is_dir() and (direct / "subtype").is_dir():
        return direct

    for candidate in extract_root.rglob(GENE_DATA_DIRNAME):
        if candidate.is_dir() and (candidate / "subtype").is_dir():
            return candidate

    raise FileNotFoundError("Could not find extracted gene_count_data/subtype directory")


def ensure_data_ready() -> BootstrapResult:
    """Ensure required input files are present in data root.

    First run downloads missing artifacts and extracts gene_count_data. Subsequent
    runs skip if data already exists.
    """
    result = BootstrapResult()
    data_root = get_data_root()

    try:
        sources = _load_sources()
        data_root.mkdir(parents=True, exist_ok=True)

        results_path = data_root / RESULTS_FILENAME
        gene_data_path = data_root / GENE_DATA_DIRNAME

        if results_path.is_file():
            result.skipped.append(RESULTS_FILENAME)
        else:
            LOGGER.info("Downloading %s to %s", RESULTS_FILENAME, results_path)
            _download_to_path(sources[RESULTS_FILENAME], results_path)
            result.downloaded.append(RESULTS_FILENAME)

        if _gene_data_ready(data_root):
            result.skipped.append(GENE_DATA_DIRNAME)
        else:
            with tempfile.TemporaryDirectory(prefix="gene_data_bootstrap_", dir=data_root) as tmp:
                tmp_root = Path(tmp)
                archive_path = tmp_root / GENE_DATA_ARCHIVE
                extract_root = tmp_root / "extract"
                extract_root.mkdir(parents=True, exist_ok=True)

                LOGGER.info("Downloading %s to %s", GENE_DATA_ARCHIVE, archive_path)
                _download_to_path(sources[GENE_DATA_ARCHIVE], archive_path)

                LOGGER.info("Extracting %s", archive_path)
                _safe_extract_tar(archive_path, extract_root)

                extracted_gene_data = _find_extracted_gene_data_dir(extract_root)
                if gene_data_path.exists():
                    shutil.rmtree(gene_data_path)
                shutil.move(str(extracted_gene_data), str(gene_data_path))
                result.downloaded.append(GENE_DATA_DIRNAME)

        result.details.append(f"Data root: {data_root}")
        return result

    except (FileNotFoundError, URLError, OSError, ValueError, tarfile.TarError) as exc:
        LOGGER.exception("Failed during data bootstrap")
        result.error = str(exc)
        result.details.append(f"Data root: {data_root}")
        return result
