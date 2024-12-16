rule fastp:
    input:
        r1="reads/pe/{sample}_R1_001.fastq.gz ",
        r2="reads/pe/{sample}_R2_001.fastq.gz "
    output:
        trimmed1="trimmed/{sample}_P_R1.fastq.gz",
        trimmed2="trimmed/{sample}_P_R2.fastq.gz",
        # Unpaired reads separately
        unpaired1="trimmed/failed/{sample}_U_R1.fastq.gz",
        unpaired2="trimmed/failed/{sample}_U_R2.fastq.gz",
        failed="trimmed/failed/{sample}.failed.fastq",
        html="fastp_report/{sample}.html",
        json="fastp_report/{sample}.json"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra="--merge"
    envmodules:
        "fastp/0.23.2-GCC-11.3.0"
    resources:
        mem_mb=12000,
        cpus_per_task=2
    shell:
        "fastp -w 16 --dont_overwrite --in1 {input.r1} --in2 {input.r2} --out1 {output.trimmed1} --out2 {output.trimmed2} --unpaired1 {output.unpaired1} --unpaired2 {output.unpaired2} --failed_out {output.failed} -j {output.json} -h {output.html}"
