#!/usr/bin/env python

import sys
import os
import subprocess
from loguru import logger
from Bio import SeqIO
import argparse


def classify_mobrecon(input_fasta, input_bed, mobrecondir, outputdir):
    if not os.path.exists(os.path.join(outputdir)):
        os.mkdir(os.path.join(outputdir))
    contigreport = os.path.join(mobrecondir, "contig_report.txt")
    mobtyper = os.path.join(mobrecondir, "mobtyper_results.txt")
    output_bed = os.path.join(
        outputdir,
        "input-mobrecon_out-intersect.sorted.bed",
    )
    with open(f"{output_bed}", "w") as f:
        pass
    if not os.path.exists(mobtyper):
        logger.info(
            "No plasmids identified by mob_recon from the provided assembly sequence"
        )
        return output_bed
    else:
        for contig in SeqIO.parse(input_fasta, "fasta"):
            # print(contig.description)
            plasmidid = (subprocess.run(
                [
                    f"grep \"{contig.description}\" {contigreport} \
                    | awk -F'\\t' '$2 == \"plasmid\" {{print}}'\
                    | cut -f3"
                ],
                capture_output=True,
                shell=True,
            ).stdout.decode().strip())
            # print(plasmidid)
            if plasmidid:
                output = subprocess.run(
                    [
                        f'grep "{contig.id}" {input_bed} \
                        | sed "s/$/\\tplasmid-contig:{plasmidid}|"$(grep {plasmidid} {mobtyper} | cut -f14)"/"'
                    ],
                    capture_output=True,
                    shell=True,
                )
                # print(output)
                if output.returncode != 0:
                    logger.error(
                        f"Error in classifying mobrecon's results! grep {contig.id} & {plasmidid}"
                    )
                    logger.error(output)
                    logger.error(output.stdout.decode())
                    logger.error(output.stderr.decode())
                    sys.exit(2)
                else:
                    with open(f"{output_bed}", "a") as f:
                        f.write(str(output.stdout.decode()))
    logger.success("Completed classifying mobrecon's output")
    return output_bed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        help="[REQUIRED] Path to mobrecon's output",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--fasta",
        "-f",
        help="[REQUIRED] Path to fasta file being analyzed",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--bed",
        "-b",
        help="[REQUIRED] Path to the input bed file to query for putative\
 mobility",
        type=os.path.abspath,
        required=True,
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Directory for the output. Default is ./maybemobile_out/",
        type=os.path.abspath,
        default="./maybemobile_out/",
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
    args = parser.parse_args()
    mobrecon_outbed = classify_mobrecon(args.fasta, args.bed, args.input,
                                        args.output)


if __name__ == "__main__":
    main()
