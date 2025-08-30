#!/usr/bin/env python

import argparse 
import logging
import os
# import requests
from tqdm import tqdm
from pathlib import Path
from typing import Iterator

from raichu.errors import NoModulesException
from raichu.retromol import read_cluster_files
from raichu.module import NRPSModule
from raichu.domain.domain import RecognitionDomain, RecognitionDomainType


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("index_a_domains")


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", required=True, help="Input directory.")
    parser.add_argument("-w", "--workdir", required=True, help="Working directory for script.")
    return parser.parse_args()


def check_existence_file(filepath: str) -> bool:
    return os.path.isfile(filepath)


def iter_files_with_extension(root_dir: str, extension: str) -> Iterator[Path]:
    """
    Recursively yield all files under `root_dir` (including subdirectories)
    that end with the given `extension`.

    :param root_dir: Directory to search in.
    :param extension: File extension to match (e.g., ".txt", ".gbk").
                      Case-insensitive, must include the dot.
    :yield: pathlib.Path objects for matching files.
    """
    root_dir = Path(root_dir)
    extension = extension.lower()
    for path in root_dir.rglob("*"):
        if path.is_file() and path.suffix.lower() == extension:
            yield path


def download_file_with_progress(download_link: str, destination: str) -> None:
    """
    Download a file from a URL with a progress bar.

    :param download_link: The URL of the file to download.
    :param destination: Local path (string) where the file will be saved.
    """
    with requests.get(download_link, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get("content-length", 0))
        block_size = 8192  # 8 KB chunks

        with open(destination, "wb") as f, tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=destination,
            ascii=True,
        ) as bar:
            for chunk in r.iter_content(chunk_size=block_size):
                if chunk:  # filter out keep-alive chunks
                    f.write(chunk)
                    bar.update(len(chunk))


def main() -> None:
    """Entry point of script."""
    args = cli()
    os.makedirs(args.workdir, exist_ok=True)

    # add file handler log file, delete first if it already exists
    log_filepath = os.path.join(args.workdir, "index_a_domains.log")
    if os.path.isfile(log_filepath):
        os.remove(log_filepath)
    handler = logging.FileHandler(log_filepath)
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)

    # download Parasect model if it doesn't exist
    # model_filepath = os.path.join(args.workdir, "parasect.pt")
    # if not check_existence_file(model_filepath):
    #     logger.info(f"Downloading Parasect model to {model_filepath}")
    #     download_link_parasect = "https://zenodo.org/records/13165500/files/model.parasect?download=1"
    #     download_file_with_progress(download_link_parasect, model_filepath)
    # else:
    #     logger.info(f"Parasect model already exists at {model_filepath}")

    # open file path and make sure we can append to it, delete file if it already exists
    out_filepath = os.path.join(args.workdir, "a_domains.fasta")
    if check_existence_file(out_filepath):
        os.remove(out_filepath)
    out_file_open = open(out_filepath, "a")

    processed_files = 0
    succeeded, failed = 0, 0
    for filepath in iter_files_with_extension(args.inputs, ".txt"):
        filename = os.path.basename(filepath)
        filename = os.path.splitext(filename)[0]
        try:
            cluster, cluster_info = read_cluster_files(filepath)
            for module in cluster.modules:
                if isinstance(module, NRPSModule):
                    for domain in module.domains:
                        if isinstance(domain, RecognitionDomain) and isinstance(domain.type, RecognitionDomainType):
                            aa_seq = domain.translation
                            mod_nr = module.id
                            header = f">{filename}|{mod_nr}|{domain.type}"
                            out_file_open.write(header + "\n")
                            out_file_open.write(aa_seq + "\n")
        except NoModulesException:
            logger.warning(f"No biosynthetic modules found in {filepath}")
            failed += 1
            continue
        except Exception as e:
            logger.error(f"Error processing {filepath}: {e}")
            failed += 1
            continue
        processed_files += 1
        succeeded += 1
        
    out_file_open.close()
  
    # Log the summary of processing results
    logger.info(f"Processing complete: {succeeded} succeeded, {failed} failed.")
    logger.info(f"A-domain sequences written to {out_filepath}")


if __name__ == "__main__":
    main()