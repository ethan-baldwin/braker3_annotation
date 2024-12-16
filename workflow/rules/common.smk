import pandas as pd

# samples = pd.read_table(config["samples"],header=0,delim_whitespace=True).set_index("sample_name", drop=False)

samples = (
    pd.read_table(config["samples"],header=0,delim_whitespace=True, dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)
