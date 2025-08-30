#!/usr/bin/env python3
import argparse
import fnmatch
import logging
import os
import shutil
import tempfile
import signal
from pathlib import Path
from zipfile import ZipFile

from tqdm import tqdm

# Import your parser directly
from raichu.retromol import gbk_to_cluster_files
from raichu.errors import NoModulesException

# GenBank filename patterns to look for
GENBANK_PATTERNS = ["*.gbk", "*.gb", "*.gbff", "*.genbank"]


# -------------------------
# Timeout infra (Unix/mac)
# -------------------------

class FunctionTimeoutError(Exception):
    """Raised when a function exceeds the allowed runtime."""
    pass


def _timeout_handler(signum, frame):
    raise FunctionTimeoutError("parse timed out")


def run_with_timeout(seconds, func, *args, **kwargs):
    """
    Run `func(*args, **kwargs)` with a SIGALRM timeout.
    Works on Unix/macOS; not supported on Windows.
    """
    original = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(seconds)
    try:
        return func(*args, **kwargs)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, original)


def find_genbank_files(root: Path):
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            for pat in GENBANK_PATTERNS:
                if fnmatch.fnmatch(name, pat):
                    yield Path(dirpath) / name
                    break


def mirror_txt_path(input_path: Path, extracted_root: Path, output_batch_root: Path) -> Path:
    """Map an extracted input file to its output .txt path under the mirrored output batch root."""
    rel = input_path.relative_to(extracted_root)
    return (output_batch_root / rel).with_suffix(".txt")


def process_one_file(inp: Path, extracted_root: Path, output_batch_root: Path, timeout_seconds: int):
    """
    Parse a single GenBank file to its mirrored .txt path.
    Returns (inp, ok, err_msg)
    """
    out_path = mirror_txt_path(inp, extracted_root, output_batch_root)
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        run_with_timeout(timeout_seconds, gbk_to_cluster_files, str(inp), str(out_path))
        return inp, True, None

    except NoModulesException:
        return inp, False, None

    except FunctionTimeoutError as e:
        logging.error("TIMEOUT: %s (>%ss)", inp, timeout_seconds)
        return inp, False, str(e)

    except Exception as e:
        logging.error("FAILED: %s (%s)", inp, e)
        return inp, False, str(e)


def process_batch_zip(zip_file: Path, output_root: Path, timeout_seconds: int) -> tuple[int, int, int]:
    """Process one region_batch_*.zip and return (total, ok, fail)."""
    logging.info("Processing %s", zip_file.name)

    output_batch_root = output_root / zip_file.stem
    output_batch_root.mkdir(parents=True, exist_ok=True)

    tmp_dir = Path(tempfile.mkdtemp(prefix=f"{zip_file.stem}_"))
    extracted_root = tmp_dir / "extracted"
    extracted_root.mkdir(parents=True, exist_ok=True)

    try:
        with ZipFile(zip_file, "r") as zf:
            zf.extractall(extracted_root)

        gbk_files = list(find_genbank_files(extracted_root))
        total = len(gbk_files)
        if total == 0:
            logging.warning("No GenBank files found in %s", zip_file.name)
            return 0, 0, 0

        ok = fail = 0
        with tqdm(total=total, desc=zip_file.name, unit="file") as pbar:
            for f in gbk_files:
                _, success, _ = process_one_file(f, extracted_root, output_batch_root, timeout_seconds)
                ok += int(success)
                if not success:
                    fail += 1
                pbar.update(1)

        logging.info("Done %s | total:%d ok:%d fail:%d -> %s",
                     zip_file.name, total, ok, fail, output_batch_root)

        return total, ok, fail

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(
        description="Mirror region_batch_*.zip contents into a new dir tree with parsed cluster .txt files (no GBKs)."
    )
    ap.add_argument("zip_dir", type=Path, help="Directory containing region_batch_*.zip")
    ap.add_argument("--pattern", default="region_batch_*.zip",
                    help="Glob for batch zips (default: region_batch_*.zip)")
    ap.add_argument("--output-root", type=Path, required=True,
                    help="Root directory to write mirrored output batches and .txt files")
    ap.add_argument("--timeout-seconds", type=int, default=5,
                    help="Max seconds to parse a single file (Unix/mac only). Default: 5")
    ap.add_argument("--log-level", default="INFO",
                    choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = ap.parse_args()

    args.output_root.mkdir(parents=True, exist_ok=True)

    # log both to console and to a file in the output root
    log_file = args.output_root / "retromol_parse.log"
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode="w", encoding="utf-8")
        ]
    )
    logging.info("Writing logs to %s", log_file)

    if not args.zip_dir.is_dir():
        logging.error("zip_dir is not a directory: %s", args.zip_dir)
        raise SystemExit(2)

    zips = sorted(args.zip_dir.glob(args.pattern))
    if not zips:
        logging.error("No ZIPs matching '%s' in %s", args.pattern, args.zip_dir)
        raise SystemExit(1)

    overall_total = overall_ok = overall_fail = 0
    for i, z in enumerate(zips, start=1):
        print(f"\n=== Batch {i}/{len(zips)}: {z.name} ===")
        total, ok, fail = process_batch_zip(z, args.output_root, args.timeout_seconds)
        overall_total += total
        overall_ok += ok
        overall_fail += fail

    summary = f"\nALL DONE | total:{overall_total} ok:{overall_ok} fail:{overall_fail}"
    print(summary)
    logging.info(summary)
    print(f"Mirrored outputs at: {args.output_root.resolve()}")
    logging.info("Mirrored outputs at: %s", args.output_root.resolve())


if __name__ == "__main__":
    main()
