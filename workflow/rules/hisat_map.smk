rule hisat_map:
    input:
        index_done=lambda wildcards: f"{genome_to_fasta_path[wildcards.genome]}/{wildcards.genome}.index.1.ht2",
        reads1=lambda wildcards: [f"trimmed/{sample}_P_R1.fastq.gz" for sample in SAMPLES[SAMPLES["genome"] == wildcards.genome]["sample_name"]],
        reads2=lambda wildcards: [f"trimmed/{sample}_P_R2.fastq.gz" for sample in SAMPLES[SAMPLES["genome"] == wildcards.genome]["sample_name"]]
        # reads1=f"{TRIMDIR}/{{sample}}_P_R1.fastq.gz",
        # reads2=f"{TRIMDIR}/{{sample}}_P_R2.fastq.gz"
    output:
        "mapped_reads/{sample}.{genome}.bam"
    params:
        index=f"fasta/{genome}.index"
    envmodules:
        "HISAT2/2.2.1-gompi-2022a",
        "SAMtools/1.17-GCC-12.2.0"
    resources:
        mem_mb=40000,
        cpus_per_task=16,
        runtime=8000
    shell:
        "hisat2 -x {params.index} -1 {input.reads1} -2 {input.reads2} -p 16 --dta | samtools view -O BAM | samtools sort --threads 16 -o {output}"

rule merge_bams:
    input:
        expand("mapped_reads/{{sample}}.{genome}.bam", sample=SAMPLES)
    output:
        "mapped_reads/{genome}.merged.bam"
    envmodules:
        "SAMtools/1.17-GCC-12.2.0"
    resources:
        mem_mb=50000,
        cpus_per_task=24,
        runtime=4000
    shell:
        "samtools merge -@ 24 {output} {input}"
