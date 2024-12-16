SAMPLES = glob_wildcards("/scratch/eab77806/genome_assembly/reads/rna/trimmed_reads/Psitt-{sample}.trim.R1.fastq.gz").sample
genome = "S_psittacina"
bam_suffix = "Psitt"

rule all:
    input:
        directory(f"busco/{genome}.genome.out"),
        directory(f"busco/{genome}.transcriptome.out")

rule index:
    input:
        f"fasta/{genome}.fa"
    output:
        done=f"fasta/{genome}.index.1.ht2"
    params:
        index=f"fasta/{genome}.index"
    envmodules:
        "HISAT2/2.2.1-gompi-2022a"
    resources:
        mem_mb=30000,
        cpus_per_task=8
    shell:
        "hisat2-build {input} {params.index} -p 8"

rule hisat_map:
    input:
        index_done=f"fasta/{genome}.index.1.ht2",
        reads1=f"/scratch/eab77806/genome_assembly/reads/rna/trimmed_reads/{bam_suffix}-{{sample}}.trim.R1.fastq.gz",
        reads2=f"/scratch/eab77806/genome_assembly/reads/rna/trimmed_reads/{bam_suffix}-{{sample}}.trim.R2.fastq.gz"
    output:
        f"mapped_reads/{bam_suffix}-{{sample}}.{genome}.bam"
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
        expand(f"mapped_reads/{bam_suffix}-{{sample}}.{genome}.bam", sample=SAMPLES)
    output:
        f"mapped_reads/{bam_suffix}.{genome}.merged.bam"
    envmodules:
        "SAMtools/1.17-GCC-12.2.0"
    resources:
        mem_mb=50000,
        cpus_per_task=24,
        runtime=4000
    shell:
        "samtools merge -@ 24 {output} {input}"

rule Braker:
    input:
        fasta=f"fasta/{genome}.fa.masked",
        prot="/scratch/eab77806/genome_assembly/braker/Viridiplantae.fa",
        bam=f"mapped_reads/{bam_suffix}.{genome}.merged.bam"
    output:
        "braker/braker.gff3"
    envmodules:
        "BRAKER/3.0.3-foss-2022a"
    params:
        species=f"{genome}"
    resources:
        mem_mb=100000,
        cpus_per_task=32,
        runtime=2000
    shell:
        """
        braker.pl --genome={input.fasta} \
         --AUGUSTUS_CONFIG_PATH=/scratch/eab77806/genome_assembly/braker/augustus/config \
         --bam={input.bam} \
         --softmasking \
         --gff3 \
         --workingdir=braker \
         --threads=32 \
         --species={params.species} \
         --useexisting \
         --prot_seq=/home/eab77806/prot/Viridiplantae.fa \
         --TSEBRA_PATH=/apps/eb/TSEBRA/1.1.2-foss-2022a/bin/
        """

rule AGAT:
    input:
        gff="braker/braker.gff3",
        fasta=f"fasta/{genome}.fa"
    output:
        f"braker/{genome}.filter.gff3.stats.out",
        longest_iso_gff=f"braker/{genome}.longest_isoforms.gff",
        longest_iso_prot=f"braker/{genome}.longest_isoforms.faa"
    envmodules:
        "AGAT/1.1.0",
        "gffread/0.12.7-GCCcore-11.3.0"
    resources:
        mem_mb=20000,
        cpus_per_task=1
    shell:
        """
        gffread {input.gff} \
        -g {input.fasta} \
        -o braker/{genome}.nofilter.gff3 \
        --keep-genes \
        -O \
        -d braker/gff_duplication_info.nofilter \
        -x braker/splicedCDS.nofilter.fasta -y braker/peptides.nofilter.fasta \
        -W -S --sort-alpha

        gffread {input.gff} \
        -g {input.fasta} \
        -o braker/{genome}.filter.gff3 \
        --keep-genes \
        -O -V -H -B -P --adj-stop \
        -M -d braker/gff_duplication_info -K \
        --force-exons --gene2exon --t-adopt \
        -x braker/splicedCDS.filter.fasta -y braker/peptides.filter.fasta \
        -W -S --sort-alpha -l 300 -J

        # get stats for preliminary gff3 file converted by AGAT
        agat_sp_statistics.pl --gff braker/{genome}.nofilter.gff3 -o braker/{genome}.nofilter.gff3.stats.out

        # get stats for final gff3 file filtered by gffread
        agat_sp_statistics.pl --gff braker/{genome}.filter.gff3 -o braker/{genome}.filter.gff3.stats.out

        # get longest isoforms for busco analysis
        agat_sp_keep_longest_isoform.pl -gff braker/{genome}.nofilter.gff3 > {output.longest_iso_gff}
        gffread {output.longest_iso_gff} -g {input.fasta} -y {output.longest_iso_prot} >/dev/null
        """
rule busco_genome:
    input:
        f"fasta/{genome}.fa"
    output:
        directory(f"busco/{genome}.genome.out")
    envmodules:
        "BUSCO/5.5.0-foss-2022a"
    resources:
        mem_mb=100000,
        cpus_per_task=24,
        runtime=4000
    shell:
        "busco -m genome -i {input} -l embryophyta_odb10 -o {output} -c 24 -f"

rule busco_transcriptome:
    input:
        f"braker/{genome}.longest_isoforms.faa"
    output:
        directory(f"busco/{genome}.transcriptome.out")
    envmodules:
        "BUSCO/5.5.0-foss-2022a"
    resources:
        mem_mb=50000,
        cpus_per_task=12,
        runtime=2000
    shell:
        "busco -m protein -i {input} -l eudicots_odb10 -o {output} -c 12"

rule entap:
    input:
        "braker/braker.aa"
    output:
        "entap_outfiles/final_results/annotated.faa"
    envmodules:
        "EnTAP/1.0.0-foss-2022a"
    resources:
        mem_mb=50000,
        cpus_per_task=12,
        runtime=4000
    shell:
        "EnTAP --runP -i {input} -t 12 -d /apps/db2/EnTAP/EggNOG_DIAMOND_Reference/bin/eggnog_proteins.dmnd -d /apps/db2/EnTAP/EggNOG_DIAMOND_Reference/bin/uniprot_sprot.dmnd --ini ~/entap_config.ini"
