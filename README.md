# zres

A tool for analysing the Z-axis distribution of amino acids using [PPM 3.0, made available by the Orientations of Proteins in Membranes database](https://opm.phar.umich.edu/).


## Usage
Given a membrane-embedded structure file, zres will return a single-row Pandas DataFrame containing a column for each amino acid type found in the transmembrane region.
Each column entry is a list of C-alpha Z-axis coordinates in angstrom relative to the membrane center Z-value. 

```
import zres
result = zres.Zres(file).run()
```

The ```batch_query.py``` automates fetching, PPM processing, and Zres analysis given a TSV file containing Uniprot IDs in a column labelled ```'Entry'```

Enjoy!