import zres
import pandas as pd

# Load in query PDB ids
query = pd.read_csv('query.csv')
query_codes = query['PDB'].to_list()

# Collect zres results
results = pd.concat([zres.Zres(i).run() for i in query_codes])
results.index = query_codes

# Save out summary
results.to_csv('results.csv')