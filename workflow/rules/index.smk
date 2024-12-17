# import shutil
#
# shutil.copyfile(src, dst)

rule index:
    input:
        fasta=lambda wildcards: f"{genome_to_fasta_path[wildcards.genome]}",
    output:
        done=lambda wildcards: f"{genome_to_fasta_path[wildcards.genome]}.index.1.ht2"
    params:
        index=lambda wildcards: f"{genome_to_fasta_path[wildcards.genome]}.index"
    envmodules:
        "HISAT2/2.2.1-gompi-2022a"
    resources:
        mem_mb=30000,
        cpus_per_task=8
    shell:
        "hisat2-build {input} {params.index} -p 8"
