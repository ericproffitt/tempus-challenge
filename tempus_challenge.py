import requests, sys

## Eric Proffitt, c. Sep, 2022

## read relevant data from VCF file
names = []
variants = []
with open("test_vcf_data.txt") as f:
	for line in f:
		if line[0] == '#':
			continue

		fields = line.strip().split('\t')

		chrom = fields[0]
		pos = fields[1]
		ref = fields[3]
		alts = fields[4].split(',')
		TC = fields[7].split('TC=')[1].split(';')[0]
		TR = fields[7].split('TR=')[1].split(';')[0].split(',')

		for alt, tr in zip(alts, TR):
			## SNV
			if len(ref) == len(alt) == 1:
				name = f'"{chrom}:g.{pos}{ref}>{alt}"'

			## insertion
			elif len(ref) == 1:
				name = f'"{chrom}:g.{pos}_{int(pos) + 1}ins{alt[1:]}"'

			## deletion
			elif len(alt) == 1:
				name = f'"{chrom}:g.{pos}_{int(pos) + len(ref)}del{ref}"'

			## indel
			else:
				name = f'"{chrom}:g.{pos}_{int(pos) + len(ref)}delins{alt}"'

			names.append(name)
			variants.append([chrom, pos, ref, alt, TC, tr, str(100 * int(tr) / int(TC))])

names = names
variants = variants

variants = dict([(name, variant) for name, variant in zip(names, variants)])

## VEP HGVS API has an upper limit of 300 variants per query
partition = [[i,i+300] for i in range(0, len(names), 300)]
partition[-1][1] = len(names)

for p in partition:
	## request data from API
	variant_str = '[' + ','.join(names[p[0]:p[1]]) + ']'

	r = requests.post("https://grch37.rest.ensembl.org/vep/human/hgvs?variant_class=1", headers={'Content-type':'application/json', 'Accept':'application/json'}, data='{{"hgvs_notations" : {}}}'.format(variant_str))
	decoded = r.json()

	## read data from API output
	for d in decoded:
		name = d['input']
		gene = d.get('transcript_consequences', [{}])[0].get('gene_symbol', 'NA')
		variant_type = d.get('variant_class', 'NA')
		variant_effect = d.get('most_severe_consequence', 'NA')
		maf = str(d.get('colocated_variants', [{}])[-1].get('minor_allele_freq', 'NA'))

		variants['"' + name + '"'] += [gene, variant_type, variant_effect, maf]

## write annotated variants to TSV file
with open("variants.tsv", 'w') as fout:
	fout.write("chrom\tpos\tref\talt\tdepth\tvar_reads\tvar_reads_perc\tgene\ttype\teffect\tmaf\n")

	for name in names:
		fout.write('\t'.join(variants[name]) + '\n')
