intermediate_dir = config["intermediate_dir"]


rule reformat_bed:
    output:
        f"{intermediate_dir}/{{batch}}/cleaned_beds/sequence_{{seqnum}}.cleaned.bed",
    input:
        fa=lambda wildcards: get_sample_path(wildcards.batch, wildcards.seqnum),
        bed=lambda wildcards: get_bed_path(wildcards.batch, wildcards.seqnum),
    conda:
        "../envs/callmemobile.yml"
    params:
        script=Path(workflow.basedir) / "scripts/format_bed.py",
    shell:
        """
        {params.script} \\
            -i {input.bed} \\
            -o {output} \\
            -f {input.fa}
        """
