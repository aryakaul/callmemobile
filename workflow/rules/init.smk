# In init.smk

import os
from pathlib import Path


configfile: "config.yaml"


def dir_input():
    return Path(config["input_dir"])


def dir_intermediate():
    return Path(config["intermediate_dir"])


def dir_output():
    return Path(config["output_dir"])


beds = {}
inputs = {}  # batch -> list of sample IDs
sample_seqnums = {}  # (batch, sample ID) -> sequence number

batch_txts = dir_input().glob("*.txt")
for x in batch_txts:
    b = os.path.basename(x)[
        :-4
    ]  # Remove '.txt' from the filename to get the batch name
    inputs[b] = []
    with open(x, "r") as f:
        seqnum = 1
        for line in f:
            sample_path = line.strip()
            inputs[b].append(sample_path)
            sample_seqnums[(b, sample_path)] = seqnum
            seqnum += 1

bed_txts = dir_input().glob("*.beds")
for x in bed_txts:
    b = os.path.basename(x)[
        :-5
    ]  # Remove '.beds' from the filename to get the batch name
    beds[b] = []
    with open(x, "r") as f:
        seqnum = 1
        for line in f:
            bed_path = line.strip()
            beds[b].append(bed_path)
            sample_seqnums[(b, bed_path)] = seqnum
            seqnum += 1


def get_batches():
    return inputs.keys()


def get_samples(_batch):
    return inputs[_batch]


def get_beds(_batch):
    return beds[_batch]


def get_sample_path(batch, seqnum):
    return inputs[batch][int(seqnum) - 1]


def get_bed_path(batch, seqnum):
    return beds[batch][int(seqnum) - 1]


def get_seqnum(_batch, _sample_id):
    return sample_seqnums[(_batch, _sample_id)]


def fn_cleanbed(batch, seqnum):
    return f"{dir_intermediate()}/{batch}/cleaned_beds/sequence_{seqnum}.cleaned.bed"


def fn_mobprocessed(batch, seqnum, sample):
    return f"{dir_intermediate()}/{batch}/mob-suite/processed/sequence_{seqnum}/{sample}/mobrecon_out/input-mobrecon_out-intersect.sorted.bed"


def fn_plasmidfinderprocessed(batch, seqnum, sample):
    # return f"{dir_intermediate()}/{batch}/mob-suite/processed/sequence_{seqnum}/{sample}/mobrecon_out/input-mobrecon_out-intersect.sorted.bed"
    return f"{dir_intermediate()}/{batch}/PlasmidFinder/processed/sequence_{seqnum}/{sample}/plasmidfinder_out.sorted.bed"


def fn_integronfinderprocessed(batch, seqnum, sample):
    return f"{dir_intermediate()}/{batch}/IntegronFinder/processed/sequence_{seqnum}/{sample}/input-ifinder_out-intersect.sorted.bed"


def fn_mefinderprocessed(batch, seqnum, sample):
    return f"{dir_intermediate()}/{batch}/mobileelementfinder/processed/sequence_{seqnum}/{sample}/input-mge_out-intersect.sorted.bed"


def fn_phigaroprocessed(batch, seqnum, sample):
    return f"{dir_intermediate()}/{batch}/phigaro/processed/sequence_{seqnum}/{sample}/input-phigaro_out-intersect.sorted.bed"


def fn_callmemobile(batch, seqnum, sample):
    return f"{dir_output()}/{batch}/sequence_{seqnum}-{sample}-callmemobile.tsv"


def fn_integronfinder_allout():
    outputs = []
    for b in get_batches():
        for s in get_samples(b):
            seqnum = get_seqnum(b, s)
            sample_base = Path(os.path.basename(s)).stem
            outputs.append(fn_integronfinderprocessed(b, seqnum, sample_base))
    return outputs


def fn_plasmidfinder_allout():
    outputs = []
    for b in get_batches():
        for s in get_samples(b):
            seqnum = get_seqnum(b, s)
            sample_base = Path(os.path.basename(s)).stem
            # outputs.append(f"{dir_intermediate()}/{b}/PlasmidFinder/raw/sequence_{seqnum}/{sample_base}/data.json")
            outputs.append(fn_plasmidfinderprocessed(b, seqnum, sample_base))
    return outputs


def fn_mobrecon_allout():
    outputs = []
    for b in get_batches():
        for s in get_samples(b):
            seqnum = get_seqnum(b, s)
            sample_base = Path(os.path.basename(s)).stem
            # outputs.append(f"{dir_intermediate()}/{b}/mob-suite/raw/sequence_{seqnum}/{sample_base}/contig_report.txt")
            outputs.append(fn_mobprocessed(b, seqnum, sample_base))
    return outputs


def fn_mobileelementfinder_allout():
    outputs = []
    for b in get_batches():
        for s in get_samples(b):
            seqnum = get_seqnum(b, s)
            sample_base = Path(os.path.basename(s)).stem
            outputs.append(fn_mefinderprocessed(b, seqnum, sample_base))
    return outputs


def fn_phispy_allout():
    outputs = []
    for b in get_batches():
        for s in get_samples(b):
            seqnum = get_seqnum(b, s)
            sample_base = Path(os.path.basename(s)).stem
            outputs.append(
                f"{dir_intermediate()}/{b}/PhiSpy/raw/sequence_{seqnum}/{sample_base}.phispy"
            )
    return outputs


def fn_phigaro_allout():
    outputs = []
    for b in get_batches():
        for s in get_samples(b):
            seqnum = get_seqnum(b, s)
            sample_base = Path(os.path.basename(s)).stem
            outputs.append(fn_phigaroprocessed(b, seqnum, sample_base))
    return outputs


def fn_callmemobile_allout():
    outputs = []
    for b in get_batches():
        for s in get_samples(b):
            seqnum = get_seqnum(b, s)
            sample_base = Path(os.path.basename(s)).stem
            outputs.append(fn_callmemobile(b, seqnum, sample_base))
    return outputs
