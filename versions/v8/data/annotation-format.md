# Domain annotation imports
Analyze FASTA before importing. Coordinates are one-based, inclusive. HMMER and InterProScan imports require exact protein identifiers and lengths but cannot verify sequence identity; ensure they were generated from the same FASTA. HMMER independent domain E-value cutoff is 1e-5, an annotation filter, not a temperature-response threshold.

JSON requires an array with the exact full protein sequence and a source:
```json
[{"id":"your_FASTA_id","sequence":"YOUR_FULL_AMINO_ACID_SEQUENCE","annotation":{"source":"Publication or annotation tool/version","domains":[{"name":"NB-ARC","start":190,"end":420,"source":"Pfam / your run","accession":"PF00931"}],"lrrRepeats":[{"start":550,"end":565},{"start":574,"end":589}]}}]
```
Replace all examples with your actual sequence and coordinates. Optional fields are `firstLrrStart`, `firstLrrEnd`, `lrrEnd`, `whdStart`, `whdEnd`, and `nlrClass`. Repeat ranges must agree with supplied first/last coordinates. A generic LRR signature span is not an individual-repeat annotation. Missing first-repeat end leaves LRR1 containment unknown. Kinase signatures are displayed only when imported or source-supported; motif Kin-2 does not establish a protein kinase domain.
