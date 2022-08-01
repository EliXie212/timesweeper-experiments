from scipy import stats
import allel


def fet_alleles(ts_list, max_allele):
    """
    Test for significant change in allele between the start and the end of
    measurement using Fisher exact test .

    Args:
        ts_list (list[np.arr]): time-serialized list of arrays of SNP or haplotype data.
                                organized with split_arr().
        max_allele???
    Returns:
        float: fisher-exact-test p-value result.
    """

    # ts_list = allel.GenotypeArray(ts_list)
    ## Get first and last geno in the sampled data
    last_genos = allel.GenotypeArray(geno_list[-1]).count_alleles(max_allele=1)
    first_genos = allel.GenotypeArray(geno_list[0]).count_alleles(max_allele=1)

    ## Calculate the changes in geno counts
    diff_genos_count = last_genos - first_genos

    ## Get the index of the major allele
    major_allele_index = np.argmax(diff_genos_count, axis=1)

    ## Get the total count for major and minor index
    major_counts = sum(diff_ts_count[major_index])
    minor_counts = sum(np.delete(diff_ts_count, major_index,axis=0))

    ## Perform fisher's exact test
    fet_or, fet_p_val = stats.fisher_exact(np.array([abs(major_counts),
                                                     abs(minor_counts)]))
                                                     
    return fet_p_val
