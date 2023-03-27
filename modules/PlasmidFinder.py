import subprocess
import os
import sys
from loguru import logger
from Bio import SeqIO


def run_plasmidfinder(input_fasta, output_path):
    plasmidfinder_output = os.path.join(
        "/".join([output_path, "plasmidfinder_out"])
    )
    if not os.path.exists(plasmidfinder_output):
        os.makedirs(plasmidfinder_output)
        logger.info(
            f"Created {plasmidfinder_output} for output of plasmidfinder"
        )
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
        # logger.debug(output.stderr.decode())
        logger.success("Completed plasmidfinder")
    return os.path.join("/".join([plasmidfinder_output, "results_tab.tsv"]))


def bedformat_plasmidfinder(pfinderout, description_to_id):
    pfinder_outputbed = os.path.dirname(pfinderout)
    pfinder_outputbed = os.path.join(
        "/".join([pfinder_outputbed, "plasmidfinder_out.sorted.bed"])
    )
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
                print(linelist[1])
                pfinder_coords = linelist[1].split("..")
                if pfinder_id in description_to_id:
                    f.write(
                        description_to_id[pfinder_id]
                        + "\t"
                        + "\t".join(pfinder_coords)
                        + "\t"
                        + "\t".join(linelist[2:])
                        + "\n"
                    )
                else:
                    f.write(
                        pfinder_id
                        + "\t"
                        + "\t".join(pfinder_coords)
                        + "\t"
                        + "\t".join(linelist[2:])
                        + "\n"
                    )
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
    output_bed = os.path.join(
        "/".join([bed, "input-pfinder_out-intersect.sorted.bed"])
    )

    for contig in SeqIO.parse(inputfasta, "fasta"):
        contig_len = len(contig)
        contig_id = contig.id
        if contig_len > max_plasmidlen:
            logger.info(
                f"{contig_id} is {contig_len:,} bp which is greater than max \
plasmid len: {max_plasmidlen:,} bp."
            )
            logger.info(f"treating {contig_id} as a chromosome.")
            logger.info(
                f"searching for input elements within {maxdist} of plasmid elements"
            )
            output = subprocess.run(
                [
                    f"closest-features --chrom {contig_id} --dist --closest {inputbed} {bedpfinder} \
                    | grep -v '|NA|NA$' \
                    | awk -F'\\t' '{{split($7, a, \"|\");if(sqrt(a[2]^2) < {maxdist}){{print $0}}}}'"
                ],
                capture_output=True,
                shell=True,
            )
            if output.returncode != 0:
                logger.error(
                    f"Error in classifying PlasmidFinder's results! closest-features {contig_id}"
                )
                logger.error(output.stdout.decode())
                logger.error(output.stderr.decode())
                sys.exit(2)
            else:
                with open(f"{output_bed}", "w") as f:
                    f.write(str(output.stdout.decode()))
                # logger.debug(output.stderr.decode())
        else:
            logger.info(f"{contig_id} is {contig_len:,} bp")
            logger.info(
                f"treating {contig_id} as a plasmid. classifying all input elements on this contig as maybe mobile"
            )
            output = subprocess.run(
                [
                    f'grep -q "{contig_id}" "{bedpfinder}" && grep "{contig_id}" "{inputbed}" \
                    | sed "s/$/\\tplasmid-contig/"'
                ],
                capture_output=True,
                shell=True,
            )
            if output.returncode != 0:
                logger.error(
                    f"Error in classifying PlasmidFinder's results! grep {contig_id}"
                )
                logger.error(output.stdout.decode())
                logger.error(output.stderr.decode())
                sys.exit(2)
            else:
                with open(f"{output_bed}", "a") as f:
                    f.write(str(output.stdout.decode()))
    logger.success("Completed classifying PlasmidFinder output")
    return output_bed
