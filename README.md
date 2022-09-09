## Tempus Challenge

The relevant data is read from *test_vcf_data.txt*.

For each variant, the variant IDs are built in HGVS sequence variant nomenclature format and submitted in query blocks of 300 to the VEP HGVS API.

Data is retrievd from the API in JSON format and all data is written to the output file, *variants.tsv*, containing the following fields,

1. chromosome
2. position
3. reference allele (humanG1Kv37)
4. alternative allele
5. depth of sequence coverage at the site of variation
6. number of reads supporting the variant
7. percentage of reads supporting the variant versus those supporting reference reads
8. gene of the variant
9. variant class
10. variant effect
11. minor allele frequency

Note that loci with multiple alternative alleles are broken up into separate rows.

Unavailable data is denoted `NA`.

CLI run command,

```bash
python3 tempus_challenge.py
```
