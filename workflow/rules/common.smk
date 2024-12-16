import pandas as pd

SAMPLES = pd.read_table(config["samples"],delim_whitespace=True).set_index("sample_name", drop=False)

# SAMPLES = (
#     pd.read_table(config["samples"],header=0,delim_whitespace=True, dtype={"sample_name": str})
#     .set_index("sample_name", drop=False)
#     .sort_index()
# )
