import gzip

gene_id_to_interval = {}
with gzip.open("/humgen/atgu1/fs03/konradk/exac/exac_browser/bundle/gencode.gtf.gz") as f:
    for line in f:
        line = line.decode('ascii').strip('\n')
        if line.startswith('#'):
            continue
        fields = line.split("\t")
        if fields[2] != "gene":
            continue
        gene_id = line.split("gene_id \"")[1].split(".")[0]
        gene_id_to_interval[gene_id] = (fields[0].replace("chr", ""), int(fields[3]), int(fields[4]))

chr22_example_intervals = [
    ['exac homeage - pcsk9', '1', 55505221, 55530525],
    ['exac homeage - ppara', '22', 46594230	, 46631326],
]

with open("top_intervals.txt", "w") as f2:
    for gene_id, chrom, start, stop in chr22_example_intervals:
        f2.write("\t".join(map(str, [gene_id, chrom, start, stop])) + "\n")
    with open("top_genes.txt") as f1:
        for line in f1:
            gene_id = line.strip('\n')
            chrom, start, stop = gene_id_to_interval[gene_id]
            f2.write("\t".join(map(str, [gene_id, chrom, start, stop])) + "\n")
