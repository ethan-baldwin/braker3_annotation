genome_to_samples = SAMPLES.groupby("genome")["sample_name"].apply(list).to_dict()
# genome_to_fasta_path = GENOMES.set_index("genome_id")["genome_path"].to_dict()

rule Braker:
    input:
        fasta=lambda wildcards: f"{genome_to_fasta_path[wildcards.genome]}/{wildcards.genome}.fa.masked",
        rna=lambda wildcards: genome_to_samples[wildcards.genome]
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
