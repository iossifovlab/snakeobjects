def run(proj,OG):
    for n in "ABC":
        OG.add("sample", n, {"fqId":n+".fastq"})
    OG.add("vcf", "o", {}, OG["sample"])
