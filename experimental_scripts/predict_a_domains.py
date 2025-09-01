#!/usr/bin/env python3

# pip install git+https://github.com/BTheDragonMaster/parasect.git@webapp

import argparse 
import logging
import os
from typing import List

import joblib
import requests
from tqdm import tqdm
from parasect.api import Result, run_parasect
from sklearn.ensemble import RandomForestClassifier


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("predict_a_domains")


def cli() -> argparse.Namespace:
    """Command line interface."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", required=True, help="Path to input FASTA file with A-domain protein sequences.")
    parser.add_argument("-w", "--workdir", required=True, help="path to working directory.")
    return parser.parse_args()


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


def fasta_batches(path: str, batch_size=1000):
    batch = []
    header = None 
    seq_parts = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # flush previous record
                if header is not None:
                    seq = "".join(seq_parts)
                    batch.extend([header, seq])
                    if len(batch) >= batch_size * 2:
                        yield batch
                        batch = []
                header = line.strip()
                seq_parts = []
            else:
                seq_parts.append(line)

        # last record
        if header is not None:
            seq = "".join(seq_parts)
            batch.extend([header, seq])
    
    if batch:
        yield batch


def predict_with_parasect(aa_seq: str, temp_dir_path: str, model: RandomForestClassifier) -> List[Result]:
    try:
        return run_parasect(aa_seq, "fasta", temp_dir_path, model)
    except Exception as e:
        logger.error(f"Error running Parasect on sequence: {e}")
        return []


def main() -> None:
    """Entry point."""
    args = cli()
    os.makedirs(args.workdir, exist_ok=True)

    # Add file handler to logger
    log_filepath = os.path.join(args.workdir, "predict_a_domains.log")
    if os.path.exists(log_filepath):
        os.remove(log_filepath)
    file_handler = logging.FileHandler(log_filepath)
    file_handler.setLevel(logging.INFO)
    logger.addHandler(file_handler)

    logger.info(f"args: {args}")

    # Download Parasect model if it doesn't exist
    model_filepath = os.path.join(args.workdir, "parasect.pt")
    if not os.path.isfile(model_filepath):
        logger.info(f"Downloading Parasect model to {model_filepath}")
        download_link_parasect = "https://zenodo.org/records/13165500/files/model.parasect?download=1"
        download_file_with_progress(download_link_parasect, model_filepath)
    else:
        logger.info(f"Parasect model already exists at {model_filepath}")

    # load model
    model = joblib.load(model_filepath)
    logger.info("Loaded Parasect model: %s", model_filepath)
    logger.info(model)

    outfile_path = os.path.join(args.workdir, "predictions.tsv")
    outfile = open(outfile_path, "w")

    # create temp dir
    temp_dir_path = os.path.join(args.workdir, "temp")
    os.makedirs(temp_dir_path, exist_ok=True)
    header_written = False

    # loop over batches, predict, write out
    for batch_idx, batch in enumerate(fasta_batches(args.fasta, batch_size=1000)):
        logger.info(f"Processing batch {batch_idx + 1} of {len(batch)//2} sequences")
        batch_fasta = "\n".join(batch)
        preds = predict_with_parasect(batch_fasta, temp_dir_path, model)
        for pred in preds:
            domain_name = pred._domain.protein_name
            domain_preds = list(zip(pred._prediction_labels, pred._predictions))
            domain_preds.sort(key=lambda x: x[0])  # sort by label to make sure they are always in same order
            if not header_written:
                labels = [x[0] for x in domain_preds]
                header = ["domain_id"] + labels
                outfile.write("\t".join(header) + "\n")
                header_written = True
            # get pred vals and make sure always 3 decimals
            pred_vals = [f"{x[1]:.3f}" for x in domain_preds]
            outfile.write(f"{domain_name}\t" + "\t".join(pred_vals) + "\n")
        
    outfile.close()
    logger.info(f"Predictions written to {outfile}")

if __name__ == "__main__":
    main()
