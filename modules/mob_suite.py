import sys
import os
import subprocess
from loguru import logger


def run_mobrecon(ifasta):
    mobrecon_output = os.path.join("/".join([output_path, "mobrecon_out"]))
    if not os.path.exists(mobrecon_output):
        os.makedirs(mobrecon_output)
        logger.info(f"Created {mobrecon_output} for output of mob_recon")
    output = subprocess.run(
        [
            "mob_recon",
            "--infile",
            f"{ifasta}",
            "--num_threads",
            f"{threads}",
            "--outdir",
            f"{mobrecon_output}",
            "--force",
        ],
        capture_output=True,
    )
    if output.returncode != 0:
        logger.error("Error in mob_recon!")
        logger.error(output.stdout.decode())
        logger.error(output.stderr.decode())
        sys.exit(2)
    else:
        logger.debug(output.stdout.decode())
        logger.success("Completed mob_recon")
    return os.path.join("/".join([mobrecon_output, "contig_report.txt"]))


def classify_mobrecon(input_fasta, input_bed, mobrecondir):
    contigreport = os.path.join("/".join([mobrecondir, "contig_report.txt"]))
    mobtyper = os.path.join("/".join([mobrecondir, "mobtyper_results.txt"]))
    if not os.path.exists(mobtyper):
        logger.info(
            "No plasmids identified by mob_recon from the provided assembly sequence"
        )
        return ""
    else:
        bed = os.path.dirname(mobrecondir)
        output_bed = os.path.join(
            "/".join([bed, "input-mobrecon_out-intersect.sorted.bed"])
        )

        for contig in SeqIO.parse(inputfasta, "fasta"):
            contig_description = contig.description
            contig_id = contig.id
            output = subprocess.run(
                [
                    f"grep -q '{contig_description}' {mobtyper} \
                    grep 'plasmid' \
                    "
                ],
                capture_output=True,
                shell=True,
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
