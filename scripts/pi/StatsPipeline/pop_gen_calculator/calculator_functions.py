from bz2 import BZ2File as bzopen
import bz2

from os import stat
from textwrap import indent
import numpy as np
import math
import pandas as pd
# import allel


def calc_counts_single_sample(ref_freq, ref_depth):
    n_nucleotide1 = int(ref_freq * ref_depth)
    n_nucleotide2 = ref_depth - n_nucleotide1
    return n_nucleotide1, n_nucleotide2


def calc_pi_single_loci(ref_freq, ref_depth, method):
    if method == "Regular":
        return ref_freq * (1 - ref_freq) * float(2)

    elif method == "Schloissnig":
        n_nucleotide1 = int(ref_freq * ref_depth)
        n_nucleotide2 = ref_depth - n_nucleotide1
        total_count_at_loci = ref_depth
        return (
            n_nucleotide1
            / total_count_at_loci
            * n_nucleotide2
            / (total_count_at_loci - 1)
        ) + (
            n_nucleotide2
            / total_count_at_loci
            * n_nucleotide1
            / (total_count_at_loci - 1)
        )


def expected_global_pi_single_loci(sample_freq_vec, sample_depth_vec, method):
    """

  args:
    method: regular or schloissnig
  """
    ref_counts = [
        int(float(freq) * float(depth))
        for freq, depth in zip(sample_freq_vec, sample_depth_vec)
    ]
    alt_counts = [
        int(depth) - int(ref_count)
        for ref_count, depth in zip(ref_counts, sample_depth_vec)
    ]
    total_N = sum(sample_depth_vec)

    if method == "Regular":
        pi_vec = [
            calc_pi_single_loci(sample_freq, sample_depth, "Regular")
            for sample_freq, sample_depth in zip(sample_freq_vec, sample_depth_vec)
        ]

        # Weighted pi by reads for each sample
        pi_within = (
            sum([(pi * depth) for pi, depth in zip(pi_vec, sample_depth_vec)]) / total_N
        )
        # Pooled reads
        pi_between = calc_pi_single_loci(
            float(sum(ref_counts)) / total_N, total_N, "Regular"
        )

    elif method == "Schloissnig":
        pi_vec = [
            calc_pi_single_loci(sample_freq, sample_depth, "Schloissnig")
            for sample_freq, sample_depth in zip(sample_freq_vec, sample_depth_vec)
        ]
        # paired of samples version otherwise change denom or change this to np.mean
        pi_within = sum(pi_vec) / 2

        # works only for two sample case
        pi_between = (
            (ref_counts[0] / sample_depth_vec[0])
            * (alt_counts[1] / sample_depth_vec[1])
        ) + (
            (
                (alt_counts[0] / sample_depth_vec[0])
                * (ref_counts[1] / sample_depth_vec[1])
            )
        )
    return pi_within, pi_between


def calc_fst(pi_within, pi_between):
    """ can supply loci pi or mean across loci pi"""

    if pi_between == 0:
        F_st = 0
    else:
        F_st = 1 - (pi_within / pi_between)
    return F_st


# def calc_Hudson_fst(allele_counts_array1, allele_counts_array2, method, summarymode):

#     ac1 = allel.AlleleCountsArray(allele_counts_array1)
#     ac2 = allel.AlleleCountsArray(allele_counts_array2)
#     num, den = allel.hudson_fst(ac1, ac2)
#     # print("num:")
#     # print(num)

#     # print("den:")
#     # print(den)

#     if method == "per_variant":
#         # F_st = [num_i / den_i if den_i != 0 else 0 for num_i, den_i in zip(num, den)]

#         if len(den) == 1 and den[0] == 0:
#             F_st = np.nan
#         else:
#             F_st = num / den

#     else:
#         if np.sum(den) == 0:
#             F_st = np.nan
#         else:
#             F_st = np.sum(num) / np.sum(den)

#     if summarymode:
#         return list(num), list(den), F_st
#     else:
#         return F_st


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))


def compute_stat_for_sample_pair(
    depth_path,
    freq_path,
    snpinfo_path,
    gene_path,
    pangenome_path,
    blacklist_path,
    sample_id_vec,
    end_index,
    calc_type,
    statistic,
    coverage_min,
    summarymode,
):
    # print("hello")
    # For summary mode
    n_total_loci = 0
    ## For fst ratio of means:
    pi_within_sum: int = 0
    pi_between_sum: int = 0

    ## For fst mean of ratios:
    fst_ratio_sum: float = 0

    # For loci mode
    loci_index = []
    loci_name = []
    cn_vec = []
    gene_vec = []
    pi_along_loci = []
    pi_within_along_loci = []
    pi_between_along_loci = []
    fst_along_loci = []
    pi_sample1 = []
    pi_sample2 = []

    depth_loci_sample1 = []
    depth_loci_sample2 = []
    freq_loci_sample1 = []
    freq_loci_sample2 = []

    allele_counts_array1 = []
    allele_counts_array2 = []
    current_gene_id = ""
    secondary_gene_id = ""

    # Load gene info

    include_gene = True
    if include_gene:
        gene_f = bzopen(gene_path, "rb")
        gene_lookup_lines = gene_f.readlines()
        gene_lookup = []
        for line in gene_lookup_lines:
            gene_lookup.append(line.decode().rstrip().split("\t"))
        gene_lookup = np.array(gene_lookup)

        # ("gene lookup shape")
        # print(gene_lookup.shape)

        gene_sample_locations = [
            list(gene_lookup[0]).index(sample_id) for sample_id in sample_id_vec
        ]
        pangenome_lookup = pd.read_csv(
            pangenome_path, compression="gzip", sep="\t", index_col=0
        )

        blacklist_genes = pd.read_csv(blacklist_path, sep=",", header=None)

        blacklist_genes = list(blacklist_genes.iloc[:, 0])

    with bzopen(depth_path, "r") as bzdepth, bzopen(freq_path, "r") as bzfreq, bzopen(
        snpinfo_path, "r"
    ) as bzsnp:
        for i, (depth_line, freq_line, snp_line) in enumerate(
            zip(bzdepth, bzfreq, bzsnp)
        ):
            # Progress update
            if i % 10000 == 0:
                print(i)

            if i == end_index:
                break
            elif i == 0:
                colnames = depth_line.decode().rstrip().split("\t")

                # print("decoded colnames")
                # print(colnames)

                if len(set(sample_id_vec) - set(colnames)) > 0:
                    ("not all samples present")
                    return ("Nothing", "Nothing", "Nothing", "Nothing", "Nothing")
                else:
                    sample_locations = [
                        colnames.index(sample_id) for sample_id in sample_id_vec
                    ]

            elif i > 0:
                # else:
                decoded_depth = depth_line.decode().rstrip().split("\t")
                decoded_freq = freq_line.decode().rstrip().split("\t")
                decoded_snpinfo = snp_line.decode().rstrip().split("\t")
                # print("decode")

                sample_depth_vec = [
                    int(decoded_depth[sample_location])
                    for sample_location in sample_locations
                ]

                sample_freq_vec = [
                    float(decoded_freq[sample_location])
                    for sample_location in sample_locations
                ]

                ## ALL FILTERS
                # 1. Sufficient
                sufficient_coverage = [
                    sample_depth >= coverage_min for sample_depth in sample_depth_vec
                ]

                if not all(sufficient_coverage):
                    continue

                # 2.  check copy number
                # copy_numbers = []
                # if not NC
                # copy_numbers = [np.nan, np.nan]
                # include_gene_info = True
                # if include_gene_info:

                if len(decoded_snpinfo) > 6:
                    gene_id = decoded_snpinfo[6]

                    # This next chunk of code is just to prevent frequently redundant looking up oof copy numbers when the gene hasn't changed

                    if gene_id == current_gene_id:
                        pass
                    elif gene_id == secondary_gene_id:
                        gene_id = current_gene_id
                    else:
                        if gene_id not in list(gene_lookup[:, 0]):
                            # print(
                            #     "not in index of gene copy number lookup checking the pangeome"
                            # )
                            # conversion
                            if gene_id not in pangenome_lookup.index:
                                continue
                            gene_id_intersection = intersection(
                                list(pangenome_lookup.loc[gene_id][1:3]),
                                list(gene_lookup[:, 0]),
                            )
                            # print("gene_id_intersection")
                            # print(gene_id_intersection)
                            if len(gene_id_intersection) == 0:
                                continue
                            else:
                                secondary_gene_id = gene_id
                                gene_id = gene_id_intersection[0]
                    
                    if gene_id != current_gene_id:
                        current_gene_id = gene_id

                        gene_index = list(gene_lookup[:, 0]).index(gene_id)
                        # print("gene_index found")
                        # print(gene_id)

                        copy_numbers = [
                            float(gene_lookup[gene_index][gene_sample_location])
                            for gene_sample_location in gene_sample_locations
                        ]

                        # print("copy numbers")
                        # print(copy_numbers)

                    sufficient_copy_number = [
                        ((CN > 0.6) and (CN < 1.2)) for CN in copy_numbers
                    ]

                    if not all(sufficient_copy_number):
                        continue
                    else:
                        # print("Sufficient")
                        pass
                        #

                else:
                    continue
                
                # NC - non coding
                # 2.5 check if gene is in blacklist genes, and skip, if true
                if gene_id in blacklist_genes:
                    print("blacklist gene")
                    continue
                # 3. check for sufficient read support if ref_freq different
                sufficient_read_support = []
                # print("runs this code")
                for c_i in range(len(sample_freq_vec)):
                    if (sample_freq_vec[c_i] > 0) and (sample_freq_vec[c_i] < 1):
                        depth_1 = int(
                            round(sample_depth_vec[c_i] * sample_freq_vec[c_i], 0)
                        )
                        depth_2 = int(
                            round(sample_depth_vec[c_i] * (1 - sample_freq_vec[c_i]), 0)
                        )

                        # ("depth 1 and 2:  " + str(depth_1) + "," + str(depth_2))

                        if (depth_1 > 1) and (depth_2) > 1:
                            sufficient_read_support.append(True)

                        else:
                            sufficient_read_support.append(False)

                    else:
                        sufficient_read_support.append(True)
                # print(sufficient_read_support)

                if all(sufficient_read_support):
                    # print("sufficient read support")

                    pass
                    # print("copy numbers")
                    # print(copy_numbers)
                    loci_index_i = i
                    loci_name_i = decoded_depth[0]
                    cn_vec.append(np.nanmean(copy_numbers))
                    gene_vec.append(gene_id)

                else:
                    # Skip this locus

                    continue

                if statistic == "Pi":
                    depth_loci_sample1.append(sample_depth_vec[0])
                    freq_loci_sample1.append(sample_freq_vec[0])

                    if len(sample_freq_vec) > 1:
                        depth_loci_sample2.append(sample_depth_vec[1])
                        freq_loci_sample2.append(sample_freq_vec[1])

                        per_sample_ref_alt = np.array(
                            [
                                calc_counts_single_sample(sample_freq, sample_depth)
                                for sample_freq, sample_depth in zip(
                                    sample_freq_vec, sample_depth_vec
                                )
                            ]
                        )
                        ref_alt_bysamples = per_sample_ref_alt.T

                        total_ref_depth = np.sum(ref_alt_bysamples[0])
                        total_sample_depth = sum(sample_depth_vec)
                        sample_freq = total_ref_depth / total_sample_depth

                        sample_freq_vec = [sample_freq]
                        sample_depth_vec = [total_sample_depth]
                    # print("CHECK")
                    pi_vec = [
                        calc_pi_single_loci(sample_freq, sample_depth, calc_type)
                        for sample_freq, sample_depth in zip(
                            sample_freq_vec, sample_depth_vec
                        )
                    ]
                    pi_along_loci.append(pi_vec[0])

                    loci_index.append(loci_index_i)
                    loci_name.append(loci_name_i)

                elif statistic == "Fst":

                    if summarymode == False:
                        depth_loci_sample1.append(sample_depth_vec[0])
                        depth_loci_sample2.append(sample_depth_vec[1])
                        freq_loci_sample1.append(sample_freq_vec[0])
                        freq_loci_sample2.append(sample_freq_vec[1])

                    if calc_type in ["Schloissnig", "Regular"]:

                        (
                            pi_within_locus,
                            pi_between_locus,
                        ) = expected_global_pi_single_loci(
                            sample_freq_vec=sample_freq_vec,
                            sample_depth_vec=sample_depth_vec,
                            method=calc_type,
                        )
                        fst_locus = calc_fst(pi_within_locus, pi_between_locus)

                        if not math.isnan(fst_locus):
                            n_total_loci += 1
                            ## For fst ratio of means
                            pi_within_sum += pi_within_locus
                            pi_between_sum += pi_between_locus
                            ## For fst mean of ratios:
                            fst_ratio_sum += fst_locus

                        if summarymode == False:
                            loci_index.append(loci_index_i)
                            loci_name.append(loci_name_i)
                            pi_within_along_loci.append(pi_within_locus)
                            pi_between_along_loci.append(pi_between_locus)
                            fst_along_loci.append(fst_locus)

                    elif calc_type == "Hudson":

                        # need all loci stats
                        per_sample_ref_alt = [
                            calc_counts_single_sample(sample_freq, sample_depth)
                            for sample_freq, sample_depth in zip(
                                sample_freq_vec, sample_depth_vec
                            )
                        ]

                        n_ref1 = per_sample_ref_alt[0][0]
                        n_alt1 = per_sample_ref_alt[0][1]
                        n_ref2 = per_sample_ref_alt[1][0]
                        n_alt2 = per_sample_ref_alt[1][1]

                        # for summaries

                        # print("[[int(n_ref1), int(n_alt1)]]")
                        # print([[int(n_ref1), int(n_alt1)]])
                        # print("[[int(n_ref2), int(n_alt2)]]")
                        # print([[int(n_ref2), int(n_alt2)]])
                        num, den, F_st = calc_Hudson_fst(
                            [[int(n_ref1), int(n_alt1)]],
                            [[int(n_ref2), int(n_alt2)]],
                            "per_variant",
                            True,
                        )
                        ## For fst ratio of means
                        if not math.isnan(F_st):
                            n_total_loci += 1

                            pi_within_sum += num[0]
                            pi_between_sum += den[0]
                            ## For fst mean of ratios:
                            fst_ratio_sum += F_st[0]

                        if summarymode == False:
                            pi_vec = [
                                calc_pi_single_loci(
                                    sample_freq, sample_depth, "Schloissnig"
                                )
                                for sample_freq, sample_depth in zip(
                                    sample_freq_vec, sample_depth_vec
                                )
                            ]
                            # print("pi vec")
                            # print(pi_vec)
                            pi_sample1.append(pi_vec[0])
                            pi_sample2.append(pi_vec[1])

                            loci_index.append(loci_index_i)
                            loci_name.append(loci_name_i)

                            allele_counts_array1.append([int(n_ref1), int(n_alt1)])
                            allele_counts_array2.append([int(n_ref2), int(n_alt2)])

        if statistic == "Pi":

            if len(sample_id_vec) > 1:
                return (
                    loci_index,
                    loci_name,
                    gene_vec,
                    pi_along_loci,
                    depth_loci_sample1,
                    depth_loci_sample2,
                    freq_loci_sample1,
                    freq_loci_sample2,
                )

            else:
                return (
                    loci_index,
                    loci_name,
                    gene_vec,
                    pi_along_loci,
                    depth_loci_sample1,
                    freq_loci_sample1,
                )

        elif statistic == "Fst":
            if calc_type in ["Schloissnig", "Regular"]:
                if summarymode:
                    Fst_ratio_of_means = calc_fst(pi_within_sum, pi_between_sum)
                    ## For fst mean of ratios:
                    Fst_mean_of_ratios = fst_ratio_sum / n_total_loci

                    return Fst_ratio_of_means, Fst_mean_of_ratios, n_total_loci
                else:
                    mean_pi_within = np.nanmean(pi_within_along_loci)
                    mean_pi_between = np.nanmean(pi_between_along_loci)
                    fst_across_loci = calc_fst(mean_pi_within, mean_pi_between)

                    return (
                        loci_index,
                        loci_name,
                        cn_vec,
                        gene_vec,
                        pi_within_along_loci,
                        pi_between_along_loci,
                        fst_along_loci,
                        mean_pi_within,
                        mean_pi_between,
                        fst_across_loci,
                        depth_loci_sample1,
                        depth_loci_sample2,
                        freq_loci_sample1,
                        freq_loci_sample2,
                    )
            else:

                Fst_ratio_of_means = pi_within_sum / pi_between_sum
                ## For fst mean of ratios:
                Fst_mean_of_ratios = fst_ratio_sum / n_total_loci

                # print("Fst ratio of means")
                # print(Fst_ratio_of_means)
                # print("Fst mean of ratios")
                # print(Fst_mean_of_ratios)

                if summarymode:
                    return Fst_ratio_of_means, Fst_mean_of_ratios, n_total_loci
                else:

                    # print("Hudson calc")

                    fst_along_loci = calc_Hudson_fst(
                        allele_counts_array1,
                        allele_counts_array2,
                        "per_variant",
                        summarymode,
                    )
                    fst_across_loci = calc_Hudson_fst(
                        allele_counts_array1,
                        allele_counts_array2,
                        "average",
                        summarymode,
                    )

                    return (
                        loci_index,
                        loci_name,
                        cn_vec,
                        gene_vec,
                        fst_along_loci,
                        fst_across_loci,
                        depth_loci_sample1,
                        depth_loci_sample2,
                        freq_loci_sample1,
                        freq_loci_sample2,
                        pi_sample1,
                        pi_sample2,
                    )

