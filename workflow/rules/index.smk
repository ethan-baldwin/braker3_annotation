def get_genome_fasta(wildcards):
    return GENOMES.loc[wildcards.genome, "masked_genome_path"]

def get_index_outputs(wildcards):
    base = GENOMES.loc[wildcards.genome, "genome_path"]
    return [f"{base}.index.{i}.ht2" for i in range(1, 9)]

rule index:
    input:
        fasta=lambda wildcards: f"{genome_to_fasta_path[wildcards.genome]}",
    output:
        done=lambda wildcards: f"{genome_to_fasta_path[wildcards.genome]}.index.1.ht2"
    envmodules:
        "HISAT2/2.2.1-gompi-2022a"
    resources:
        mem_mb=30000,
        cpus_per_task=8
    shell:
        "hisat2-build {input} -p 8"
