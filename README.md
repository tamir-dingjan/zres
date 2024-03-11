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


The ```aligned_analysis.py``` allows using a multiple sequence alignment to specify a threshold position in the alignment after which to analyse residues for their Z-coordinate. To run the analysis beginning from, for example, position 46 in the sequence alignment, run the following:
```python aligned_analysis.py -a alignment.fasta -s 46```

Note that the sequence alignment must be in ```.fasta``` formatting, with sequence IDs located following a prefix containing ```>.*|```, e.g.,
```>tr|ID1234|following_text_not_considered```

Enjoy!