
rule fastp:
    input:
        r1=f"{READSDIR}/{{sample}}_R1_001.fastq.gz",
        r2=f"{READSDIR}/{{sample}}_R2_001.fastq.gz"
    output:
        trimmed1=f"{TRIMDIR}/{{sample}}_P_R1.fastq.gz",
        trimmed2=f"{TRIMDIR}/{{sample}}_P_R2.fastq.gz",
        # Unpaired reads separately
        unpaired1=f"{TRIMDIR}/failed/{{sample}}_U_R1.fastq.gz",
        unpaired2=f"{TRIMDIR}/failed/{{sample}}_U_R2.fastq.gz",
        failed=f"{TRIMDIR}/failed/{{sample}}.failed.fastq",
        html="fastp_report/{sample}.html",
        json="fastp_report/{sample}.json"
    envmodules:
        "fastp/0.23.2-GCC-11.3.0"
    resources:
        mem_mb=12000,
        cpus_per_task=16
    shell:
        "fastp -w {resources.cpus_per_task} --dont_overwrite --in1 {input.r1} --in2 {input.r2} --out1 {output.trimmed1} --out2 {output.trimmed2} --unpaired1 {output.unpaired1} --unpaired2 {output.unpaired2} --failed_out {output.failed} -j {output.json} -h {output.html}"
