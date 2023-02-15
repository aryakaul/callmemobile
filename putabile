#!/usr/bin/env python3

import argparse
import subprocess
import multiprocessing
import os
import glob
from loguru import logger
import sys
import pandas as pd


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
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        logger.debug(output.stdout.decode())
        logger.debug(output.stderr.decode())
        logger.success("Completed mobileelementfinder")
    return "".join([mefinder_output, ".csv"])


def run_integronfinder(input_fasta):
    integronfinder_output = os.path.join("/".join([output_path, "integronfinder_out"]))
    if not os.path.exists(integronfinder_output):
        os.makedirs(integronfinder_output)
        logger.info(f"Created {integronfinder_output} for output of Integron-Finder")
    output = subprocess.run(
        [
            "integron_finder",
            "--local-max",
            "--cpu",
            f"{threads}",
            "--circ",
            "--outdir",
            f"{integronfinder_output}",
            f"{input_fasta}",
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in IntegronFinder!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        logger.debug(output.stdout.decode())
        logger.debug(output.stderr.decode())
        logger.success("Completed IntegronFinder")
    return glob.glob(os.path.join("/".join([output_path, "integronfinder_out", "*.integrons"])))[0]

def bedoutput_integronfinder(ifinder_out):
    return ""

def run_plasmidfinder(input_fasta):
    plasmidfinder_output = os.path.join("/".join([output_path, "plasmidfinder_out"]))
    if not os.path.exists(plasmidfinder_output):
        os.makedirs(plasmidfinder_output)
        logger.info(f"Created {plasmidfinder_output} for output of plasmidfinder")
    output = subprocess.run(
        [
            "plasmidfinder.py",
            "-i",
            f"{input_fasta}",
            "-o",
            f"{plasmidfinder_output}",
            "-x",
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in plasmidfinder!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        logger.debug(output.stdout.decode())
        logger.debug(output.stderr.decode())
        logger.success("Completed plasmidfinder")


def bedformat_mobileelementfinder(mge_outputcsv):
    mge_outputbed = os.path.dirname(mge_outputcsv)
    mge_outputbed = os.path.join("/".join([mge_outputbed, "mge_out.sorted.bed"]))
    output = subprocess.run(
        [
            f"grep -v '^#' {mge_outputcsv} \
            | csvtk cut -f contig,start,end,name,type -T \
            | sed 1d \
            | awk -F'\\t' '{{split($1, a, \" \"); printf \"%s\t\", a[1]; for(i=2;i<=NF;++i) if(i!=NF){{printf \"%s\\t\", $i}}else{{printf \"%s\\n\", $i}}}}' \
            | sort-bed -",
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error("Error in formatting mge output as bedfile!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{mge_outputbed}", "w") as f:
            f.write(str(output.stdout.decode()))
        logger.debug(output.stderr.decode())
        logger.success("Completed reformatting mge output")
        return mge_outputbed


def classify_mobileelementfinder(inputbed, bedmge):
    maxdist = 20000

    bed = os.path.dirname(bedmge)
    output_bed = os.path.join("/".join([bed, "input-mge_out-intersect.sorted.bed"]))
    # output = subprocess.run(
    # [
    # f"bedops --element-of 1 {inputbed} {bedmge}"
    # ],
    # capture_output=True,
    # shell=True,
    # )
    # if output.returncode != 0:
    # logger.error("Error in classifying mge's results! bedops --element-of")
    # logger.error(output.stdout.decode())
    # logger.error(output.stderr.decode())
    # sys.exit(2)
    # else:
    # with open(f"{output_bed}", "w") as f:
    # f.write(str(output.stdout.decode()))
    # logger.debug(output.stderr.decode())

    output = subprocess.run(
        [
            f"closest-features --dist --closest {inputbed} {bedmge} \
            | awk -F'\\t' '{{split($7, a, \"|\"); if(a[2] < {maxdist}) {{print $0}}}}'"
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error("Error in classifying mge's results! closest-features")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{output_bed}", "w") as f:
            f.write(str(output.stdout.decode()))
        logger.debug(output.stderr.decode())
    logger.success("Completed classifying mge output")
    return output_bed


def run_classify_mobileelementfinder(input_fasta, input_bed):
    mgeout = run_mobileelementfinder(input_fasta)
    bedmgeout = bedformat_mobileelementfinder(mgeout)
    classify_mobileelementfinder(input_bed, bedmgeout)

def run_classify_integronfinder(input_fasta, input_bed):
    ifinderout = run_integronfinder(input_fasta)
    bedifinder = bedformat_integronfinder(ifinderout)
    classify_integronfinder(input_bed, bedifinder)


def classify_regions_plasmidfinder():
    return "null"

def run_tools(input_fasta, input_bed):

    # run mobileelementfinder
    run_classify_mobileelementfinder(input_fasta, input_bed)

    # run integron_finder
    run_classify_integronfinder(input_fasta, input_bed)

    # run plasmidfinder.py
    run_plasmidfinder(input_fasta)


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

    # run all tools
    run_tools(args.input, args.bed)


if __name__ == "__main__":
    main()
