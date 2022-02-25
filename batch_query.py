import zres
import pandas as pd
import numpy as np

# Load in query PDB ids
query = pd.read_csv('query.csv')
query_codes = query['PDB'].to_list()

valid_results = []
valid_codes = []
# Collect zres results
for i in query_codes:
    analysis = zres.Zres(i).run()
    if not analysis is None:
        valid_results.append(analysis)
        valid_codes.append(i)

results = pd.concat(valid_results)
results.index = valid_codes

# Save out summary
results.to_csv('results.csv')