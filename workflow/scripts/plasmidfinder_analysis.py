#!/usr/bin/env python
import argparse
import subprocess
import os
import sys
from loguru import logger
from Bio import SeqIO


def check_input(input_fasta, input_bed):
    fasta_sequences = SeqIO.parse(open(input_fasta), "fasta")
    # define a dictionary of chromosome keys and values
    description_to_id = {}
    for i in fasta_sequences:
        description_to_id[i.description] = i.id
    logger.debug(f"Fasta header to identifier is: {description_to_id}")
    output_bed = os.path.join(
        os.path.dirname(input_bed),
        os.path.basename(input_bed) + "-chromosomeid",
    )

    exit = False
    # open the input BED file for reading
    with open(input_bed, "r") as input_file:
        # open the output BED file for writing
        with open(output_bed, "w") as output_file:
            # loop through each line in the input BED file
            for line in input_file:
                # split the line into columns
                fields = line.strip().split("\t")
                if len(fields) != 4:
                    print(fields)
                    logger.error(
                        f"Number of columns is {len(fields)}! This is not 4! Please make it 4 columns!"
                    )
                    exit = True

                # check if the chromosome in the line matches a key in the dictionary
                if fields[0] in description_to_id:
                    exit = True
                    logger.error(
                        f"{fields[0]} found in bed file! Please use {description_to_id[fields[0]]} instead! Reformatted bedfile outputted here {output_bed}. Please rerun mamo with that bed file instead!"
                    )
                    # if there is a match, replace the chromosome value with the dictionary value
                    fields[0] = description_to_id[fields[0]]

                # write the updated line to the output BED file
                output_file.write("\t".join(fields) + "\n")
    if exit:
        sys.exit(2)
    else:
        os.remove(output_bed)
    return description_to_id


def bedformat_plasmidfinder(pfinderout, output, description_to_id):
    pfinder_outputbed = os.path.dirname(output)
    pfinder_outputbed = os.path.join(pfinder_outputbed,
                                     "plasmidfinder_out.sorted.bed")
    output = subprocess.run(
        [
            f"csvtk cut -t -f Contig,'Position in contig',Plasmid,Note {pfinderout} \
            | sed 1d"
        ],
        capture_output=True,
        shell=True,
    )
    if output.returncode != 0:
        logger.error("Error in formatting plasmidfinder output as bedfile!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        with open(f"{pfinder_outputbed}.unsorted", "w") as f:
            for lines in output.stdout.decode().split("\n"):
                linelist = lines.split("\t")
                if len(linelist) <= 1:
                    continue
                pfinder_id = linelist[0]
                pfinder_coords = linelist[1].split("..")
                if pfinder_id in description_to_id:
                    f.write(description_to_id[pfinder_id] + "\t" +
                            "\t".join(pfinder_coords) + "\t" +
                            "\t".join(linelist[2:]) + "\n")
                else:
                    f.write(pfinder_id + "\t" + "\t".join(pfinder_coords) +
                            "\t" + "\t".join(linelist[2:]) + "\n")
        output = subprocess.run(
            [f"sort-bed {pfinder_outputbed}.unsorted > {pfinder_outputbed}"],
            check=True,
            shell=True,
        )
        logger.success(
            f"Completed reformatting plasmidfinder output to {pfinder_outputbed}"
        )
        return pfinder_outputbed


def classify_plasmidfinder(inputfasta, inputbed, bedpfinder):
    max_plasmidlen = 200000
    maxdist = 20000

    bed = os.path.dirname(bedpfinder)
    output_bed = os.path.join(bed, "input-pfinder_out-intersect.sorted.bed")
    with open(f"{output_bed}", "w") as f:
        pass

    for contig in SeqIO.parse(inputfasta, "fasta"):
        output = subprocess.run(["grep", contig.id, bedpfinder],
                                stdout=subprocess.PIPE).stdout.decode("utf-8")
        if not output:
            logger.info(
                f"No plasmid elements found on {contig.id}. Continuing")
            continue
        output = subprocess.run(["grep", contig.id, inputbed],
                                stdout=subprocess.PIPE).stdout.decode("utf-8")
        if not output:
            logger.info(f"No tested elements are on {contig.id}. Continuing")
            continue
        contig_len = len(contig)
        if contig_len > max_plasmidlen:
            logger.info(
                f"{contig.id} is {contig_len:,} bp which is greater than max \
plasmid len: {max_plasmidlen:,} bp.")
            logger.info(f"treating {contig.id} as a chromosome.")
            logger.info(
                f"searching for input elements within {maxdist} of plasmid elements"
            )
            output = subprocess.run(
                [
                    f"closest-features --chrom {contig.id} --dist --closest {inputbed} {bedpfinder} \
                    | grep -v '|NA|NA$' \
                    | awk -F'\\t' '{{split($7, a, \"|\");if(sqrt(a[2]^2) < {maxdist}){{print $0}}}}'"
                ],
                capture_output=True,
                shell=True,
            )
            if output.returncode != 0:
                logger.error(
                    f"Error in classifying PlasmidFinder's results! closest-features {contig.id}"
                )
                logger.error(output.stdout.decode())
                logger.error(output.stderr.decode())
                sys.exit(2)
            else:
                with open(f"{output_bed}", "w") as f:
                    f.write(str(output.stdout.decode()))
                # logger.debug(output.stderr.decode())
        else:
            logger.info(f"{contig.id} is {contig_len:,} bp")
            logger.info(
                f"treating {contig.id} as a plasmid. classifying all input elements on this contig as maybe mobile"
            )
            output = subprocess.run(
                [
                    # f"grep {contig.id} {bedpfinder} | cut -f4"
                    f'grep {contig.id} {inputbed} | sed "s~$~\\tplasmid-contig|$(printf \'%q\' "$(grep {contig.id} {bedpfinder} | cut -f4)")~"'
                    # | sed "s/$/\\tplasmid-contig|$(grep {contig.id} {bedpfinder} | cut -f4 | sed \'s/$/-/\')/"'
                    # | sed "s/$/\tplasmid-contig|$(grep {contig.id} {bedpfinder} | cut -f4 | tr "\\n" "-" | sed \'s/-$/{chr(92)}n/\')/"'
                ],
                capture_output=True,
                shell=True,
            )
            logger.debug(output)
            if output.returncode != 0:
                logger.error(
                    f"Error in classifying PlasmidFinder's results! grep {contig.id}"
                )
                logger.error(output)
                logger.error(output.stdout.decode())
                logger.error(output.stderr.decode())
                sys.exit(2)
            else:
                with open(f"{output_bed}", "a") as f:
                    f.write(str(output.stdout.decode()))
    logger.success("Completed classifying PlasmidFinder output")
    return output_bed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        "-i",
        help="[REQUIRED] Path to plasmidfinder's output",
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
    args = parser.parse_args()

    description_to_id = check_input(args.fasta, args.bed)
    plasmidfinder_outbed = bedformat_plasmidfinder(args.input, args.output,
                                                   description_to_id)


if __name__ == "__main__":
    main()
