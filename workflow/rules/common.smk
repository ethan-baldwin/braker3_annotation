import pandas as pd

READSDIR=config['raw_reads_directory']
TRIMDIR = "trimmed"

SAMPLES = pd.read_table(config["samples"],header=0,delim_whitespace=True).set_index("sample_name", drop=False)
GENOMES = pd.read_table(config["genomes"],header=0,delim_whitespace=True).set_index("genome_id", drop=False)

genome_to_masked_path = GENOMES.set_index("genome_id")["masked_genome_path"].to_dict()
genome_to_fasta_path = GENOMES.set_index("genome_id")["fasta_path"].to_dict()

# genome_to_samples = SAMPLES.groupby("genome")["sample_name"].apply(list).to_dict()

# genome_to_fasta = GENOMES.set_index("genome_id")["masked_genome_path"].to_dict()
