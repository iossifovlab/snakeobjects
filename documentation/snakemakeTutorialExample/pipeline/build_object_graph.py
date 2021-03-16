def run(proj,OG):
    for n in "ABC":
        OG.add("sample", n)
    OG.add("vcf", "o", {}, OG["sample"])
