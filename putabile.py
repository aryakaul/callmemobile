#!/usr/bin/env python3

import argparse
import subprocess
import multiprocessing
import os
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
        help="[REQUIRED] Path to the input bed file denoting regions to query for putative mobility",
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

def run_mobileelementfinder():
    mefinder_output=os.path.join([output_path, 'mobileElementFinder_out'])


    # "mefinder find --contig {input_fasta} --threads {threads} {output_name}"

# def run_integronfinder():
    # print("hi")


def main():
    # set global vars
    global pipeline
    global repo_path
    global output_path
    pipeline = "PUTABILE"

    # set paths
    repo_path = os.path.dirname(os.path.realpath(__file__))

    # parse args
    args = parse_arguments(get_version())

    output_path = args.output
    threads = args.threads

    # run mobileelementfinder
    # run_mobileelementfinder()

    # run integron_finder 
    # run_integronfinder()

    # run plasmidfinder.py
    # run_plasmidfinder()


if __name__ == "__main__":
    main()
