import pandas as pd

samples = (
    pd.read_table(config["samples"], dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)
