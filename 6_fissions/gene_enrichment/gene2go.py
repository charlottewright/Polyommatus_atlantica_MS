import sys

with open(sys.argv[1], 'r') as tsv:
	gene2go = {}
	for line in tsv:
		cols = line.rstrip("\n").split("\t")
		gene = cols[0]
		go_list = []
		for col in cols:
			if col.startswith("GO"):
				go_list += col.split("|")
		try:
			gene2go[gene] += list(set(go_list))
		except KeyError:
			gene2go[gene] = list(set(go_list))

for gene, go_list in gene2go.items():
	print(gene + "\t" + ",".join(set(go_list)))
