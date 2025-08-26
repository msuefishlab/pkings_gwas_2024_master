import allel

# 1. Load genotypes
callset = allel.read_vcf('BP1.thinned.vcf.gz', fields=['calldata/GT'])
gt = allel.GenotypeArray(callset['calldata/GT'])

# 2. Count alleles per site → AlleleCountsArray with shape (n_sites, 2)
ac = gt.count_alleles()  

# 3. Folded SFS on the 2D counts array
sfs = allel.sfs_folded(ac)

print("Folded SFS bins:", sfs)

with open('BP1.folded.sfs', 'w') as out:
    out.write(' '.join(map(str, sfs)) + '\n')