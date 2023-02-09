#!/usr/bin/env python3

import argparse
import subprocess
import multiprocessing
import os
from loguru import logger
import sys
import pandas as pd
from yaml import dump as yaml_dump


def get_version():
    version = open(os.path.join(repo_path, "VERSION"), "r").readline().strip()
    return version


def parse_arguments(version):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        help="[REQUIRED] Path to the input fasta reference file.",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--bed",
        "-b",
        help="[REQUIRED] Path to the input bed file to query for putative mobility",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Directory for the output. Default is ./putativemobile_out/",
        type=os.path.abspath,
        default="./putativemobile_out/",
        required=False,
    )
    parser.add_argument(
        "--threads",
        "-t",
        help="Number of threads to use. Default is all available.",
        type=int,
        default=len(os.sched_getaffinity(0)),
        required=False,
    )
    parser.add_argument(
        "--version",
        "-V",
        help="Print program version",
        action="version",
        version=f"putabile {version}",
        # version="hi",
    )
    return parser.parse_args()


def run_mobileelementfinder(input_fasta):
    mefinder_output = os.path.join("/".join([output_path, "mobileElementFinder_out"]))
    if not os.path.exists(mefinder_output):
        os.makedirs(mefinder_output)
        logger.info(f"Created {mefinder_output} for output of mobileElementFinder")
    mefinder_output = os.path.join("/".join([mefinder_output, "out"]))
    output = subprocess.run(
        [
            "mefinder",
            "find",
            "--contig",
            f"{input_fasta}",
            "--threads",
            f"{threads}",
            f"{mefinder_output}",
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in mobileelementfinder!")
        logger.error(output.stdout)
        logger.error(output.stderr)
        sys.exit(2)
    else:
        logger.debug(output.stdout)
        logger.debug(output.stderr)
        logger.success("Completed mobileelementfinder")


def run_integronfinder(input_fasta):
    integronfinder_output = os.path.join("/".join([output_path, "integronfinder_out"]))
    if not os.path.exists(integronfinder_output):
        os.makedirs(integronfinder_output)
        logger.info(f"Created {integronfinder_output} for output of Integron-Finder")
    # mefinder_output = os.path.join("/".join([mefinder_output, "out"]))
    output = subprocess.run(
        [
            "integron_finder",
            "--local-max",
            "--cpu",
            f"{threads}",
            "--circ",
            "--outdir",
            f"{integronfinder_output}",
            f"{input_fasta}"
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in IntegronFinder!")
        logger.error(output.stdout)
        logger.error(output.stderr)
        sys.exit(2)
    else:
        logger.debug(output.stdout)
        logger.debug(output.stderr)
        logger.success("Completed IntegronFinder")


def main():
    # set global vars
    global pipeline
    global repo_path
    global output_path
    global threads
    pipeline = "PUTABILE"

    # set paths
    repo_path = os.path.dirname(os.path.realpath(__file__))

    # parse args
    args = parse_arguments(get_version())

    output_path = args.output
    threads = args.threads

    # run mobileelementfinder
    run_mobileelementfinder(args.input)

    # run integron_finder
    run_integronfinder(args.input)

    # run plasmidfinder.py
    # run_plasmidfinder()


if __name__ == "__main__":
    main()
