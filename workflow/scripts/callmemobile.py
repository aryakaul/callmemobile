#!/usr/bin/env python

import argparse
import sys
from loguru import logger


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Aggregate mobile element predictions from various tools."
    )
    parser.add_argument(
        "--integronfinder",
        help="Path to IntegronFinder output BED file",
        required=False,
    )
    parser.add_argument(
        "--plasmidfinder",
        help="Path to PlasmidFinder output BED file",
        required=False,
    )
    parser.add_argument(
        "--mob_suite", help="Path to mob-suite output BED file", required=False
    )
    parser.add_argument(
        "--phigaro", help="Path to phigaro output BED file", required=False
    )
    parser.add_argument(
        "--mobileelementfinder",
        help="Path to mobileelementfinder output BED file",
        required=False,
    )
    parser.add_argument(
        "--bed",
        help="[REQUIRED] Path to input BED file of regions of interest",
        required=True,
    )
    parser.add_argument(
        "--output", "-o", help="[REQUIRED] Path to output file", required=True
    )
    return parser.parse_args()


def read_bed_file(bed_file_path):
    regions = []
    with open(bed_file_path, "r") as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split("\t")
            if len(fields) < 3:
                logger.error(f"Invalid BED line in {bed_file_path}: {line.strip()}")
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else ""
            other_fields = fields[4:] if len(fields) > 4 else []
            regions.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "name": name,
                    "other_fields": other_fields,
                }
            )
    return regions


def regions_overlap(region1, region2):
    if region1["chrom"] != region2["chrom"]:
        return False
    # Check if regions overlap
    return not (region1["end"] <= region2["start"] or region2["end"] <= region1["start"])


def main():
    args = parse_arguments()

    # Check that at least one tool output is provided
    if not any(
        [
            args.integronfinder,
            args.plasmidfinder,
            args.mob_suite,
            args.phigaro,
            args.mobileelementfinder,
        ]
    ):
        logger.error("At least one tool output must be provided.")
        sys.exit(1)

    # Read input BED file
    logger.info(f"Reading input BED file: {args.bed}")
    input_regions = read_bed_file(args.bed)

    if not input_regions:
        logger.error("No valid regions found in input BED file.")
        sys.exit(1)

    # Initialize tool outputs
    tool_bed_regions = {}
    if args.integronfinder:
        logger.info(f"Reading IntegronFinder BED file: {args.integronfinder}")
        tool_bed_regions["IntegronFinder"] = read_bed_file(args.integronfinder)
    if args.plasmidfinder:
        logger.info(f"Reading PlasmidFinder BED file: {args.plasmidfinder}")
        tool_bed_regions["PlasmidFinder"] = read_bed_file(args.plasmidfinder)
    if args.mob_suite:
        logger.info(f"Reading mob-suite BED file: {args.mob_suite}")
        tool_bed_regions["mob_suite"] = read_bed_file(args.mob_suite)
    if args.phigaro:
        logger.info(f"Reading phigaro BED file: {args.phigaro}")
        tool_bed_regions["phigaro"] = read_bed_file(args.phigaro)
    if args.mobileelementfinder:
        logger.info(
            f"Reading mobileelementfinder BED file: {args.mobileelementfinder}"
        )
        tool_bed_regions["mobileelementfinder"] = read_bed_file(
            args.mobileelementfinder
        )

    tools = list(tool_bed_regions.keys())

    # Initialize overlap dictionary
    overlap_dict = {}
    for region in input_regions:
        key = (region["chrom"], region["start"], region["end"], region["name"])
        overlap_dict[key] = {tool_name: None for tool_name in tools}

    # For each tool, perform overlap check and populate overlap_dict
    for tool_name, tool_regions in tool_bed_regions.items():
        logger.info(f"Processing overlaps with {tool_name}")
        for input_region in input_regions:
            key = (
                input_region["chrom"],
                input_region["start"],
                input_region["end"],
                input_region["name"],
            )
            for tool_region in tool_regions:
                if regions_overlap(input_region, tool_region):
                    # Get annotation from tool region
                    annotation = tool_region["name"]
                    if tool_region["other_fields"]:
                        annotation += "|" + "|".join(tool_region["other_fields"])
                    if overlap_dict[key][tool_name]:
                        overlap_dict[key][tool_name] += f";{annotation}"
                    else:
                        overlap_dict[key][tool_name] = annotation
                    # Break if overlap is found to prevent multiple counts
                    # Remove the above line if multiple overlaps per region are possible
                    # break

    # Output results
    with open(args.output, "w") as out_f:
        header = ["#chrom", "start", "end", "name"] + tools
        out_f.write("\t".join(header) + "\n")
        for region in input_regions:
            key = (
                region["chrom"],
                region["start"],
                region["end"],
                region["name"],
            )
            values = [
                region["chrom"],
                str(region["start"]),
                str(region["end"]),
                region["name"],
            ]
            for tool_name in tools:
                value = overlap_dict[key][tool_name] if overlap_dict[key][tool_name] else "no"
                values.append(value)
            out_f.write("\t".join(values) + "\n")
    logger.info(f"Aggregation complete. Output written to {args.output}")


if __name__ == "__main__":
    main()


