# zres

A tool for analysing the Z-axis distribution of amino acids in the OPM database.


## Usage
Given a PDB code, zres will return a single-row Pandas DataFrame containing a column for each amino acid type found in the transmembrane region.
Each column entry is a list of C-alpha Z-axis coordinates normalised between 0 (for the inner membrane barrier) and 1 (for the outer membrane barrier).

```
import zres
result = zres.Zres(PDB_CODE).run()
```
