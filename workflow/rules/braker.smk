genome_to_samples = SAMPLES.groupby("genome")["sample_name"].apply(list).to_dict()
genome_to_fasta = GENOMES.set_index("genome_id")["masked_genome_path"].to_dict()

rule Braker:
    input:
        fasta=lambda wildcards: genome_to_fasta[wildcards.genome],
        rna=lambda wildcards: genome_to_samples[wildcards.genome],
        trimmed1=f"{TRIMDIR}/{{sample}}_P_R1.fastq.gz",
        trimmed2=f"{TRIMDIR}/{{sample}}_P_R2.fastq.gz"
    output:
        "{genome}/braker/braker.gff3"
    envmodules:
        "BRAKER/3.0.3-foss-2022a"
    params:
        species="{genome}",
        prot=config['prot'],
        tsebra_path=config['tsebra_path']
    resources:
        mem_mb=100000,
        cpus_per_task=32,
        runtime=2000
    shell:
        """
        braker.pl --genome={input.fasta} \
         --AUGUSTUS_CONFIG_PATH=/scratch/eab77806/genome_assembly/braker/augustus/config \
         --rnaseq_sets_ids={input.rna} \
         --rnaseq_sets_dirs={READSDIR} \
         --softmasking \
         --gff3 \
         --workingdir=braker \
         --threads=32 \
         --species={params.species} \
         --useexisting \
         --prot_seq={params.prot} \
         --TSEBRA_PATH={params.tsebra_path}
        """
