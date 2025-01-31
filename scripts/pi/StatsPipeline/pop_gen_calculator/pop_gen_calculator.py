import argparse, sys, os
from datetime import datetime
import yaml
import pandas as pd
import numpy as np
import calculator_functions
import bz2


def main():
    start = datetime.now()
    print(start)

    if local:
        config_file_path = "~/Wasney-Briscoe/scripts/pi/StatsPipeline/configs/config.yaml"
    else:

        config_file_path = "~/Wasney-Briscoe/scripts/pi/StatsPipeline/configs/config.yaml"

    with open(config_file_path, "r") as f:
        config = yaml.safe_load(f)
    # print(config)

    if local:
        snp_dir = config[study]["local_snp_dir"]
        input_dir = config[study]["local_input_dir"]
        gene_dir = config[study]["local_gene_dir"]
        pangenome_dir = config["pangenome"]["local"]
        blacklist_path = config["blacklist"]["local"]

    else:
        snp_dir = config[study]["snp_dir"]
        input_dir = config[study]["input_dir"]
        gene_dir = config[study]["gene_dir"]
        pangenome_dir = config["pangenome"]["hoffman"]
        blacklist_path = config["blacklist"]["hoffman"]

    depth_path = snp_dir + strain + "/" + "snps_depth.txt.bz2"
    freq_path = snp_dir + strain + "/" + "snps_ref_freq.txt.bz2"
    snpinfo_path = snp_dir + strain + "/" + "snps_info.txt.bz2"
    gene_path = gene_dir + strain + "/" + "genes_copynum.txt.bz2"
    pangenome_path = pangenome_dir + strain + "/" + "gene_info.txt.gz"

    print("Sample1")
    print(sample_id)

    print("Sample2")
    print(sample_id2)

#     output_dir = (                                                                        #MW 12/12/23: directing single sample pi elsewhere
#         input_dir + "Calculated_" + calc_statistic + "_MinCoverage" + str(cov_min) + "/"
#     )

    output_dir = (
        input_dir + "Single" + calc_statistic + "_MinCoverage" + str(cov_min) + "/"
    )
#     if not os.path.exists(output_dir): #MW 12/13/24: Made this change, see altered script below
#         os.makedirs(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    ## Arrangement

    print("paired pi")

    input_sample_id_vec = [sample_id, sample_id2]

    if sample_id2 == "":
        input_sample_id_vec = [sample_id]

    if calc_statistic == "Pi":

        if sample_id2 == "":

            (
                loci_index,
                loci_name,
                gene_vec,
                pi_along_loci,
                depth_loci_sample1,
                freq_loci_sample1,
            ) = calculator_functions.compute_stat_for_sample_pair(
                depth_path=depth_path,
                freq_path=freq_path,
                snpinfo_path=snpinfo_path,
                gene_path=gene_path,
                pangenome_path=pangenome_path,
                blacklist_path=blacklist_path,
                sample_id_vec=input_sample_id_vec,
                end_index=end_index,
                calc_type=calc_type,
                statistic=calc_statistic,
                coverage_min=cov_min,
                summarymode=summarymode,
            )

            if pi_along_loci[0] == "Nothing":
                exit()
            print("writing to table")
            tests_pd = pd.DataFrame(
                np.nan,
                index=list(range(len(pi_along_loci))),
                columns=[
                    "loci_index",
                    "loci_name",
                    "gene_id",
                    "pi_along_loci",
                    "depth_loci_sample1",
                    "freq_loci_sample1",
                ],
            )
            tests_pd.iloc[:, 0] = loci_index
            tests_pd.iloc[:, 1] = loci_name
            tests_pd.iloc[:, 2] = gene_vec
            tests_pd.iloc[:, 3] = pi_along_loci
            tests_pd.iloc[:, 4] = depth_loci_sample1
            tests_pd.iloc[:, 5] = freq_loci_sample1

            genomewide_pi = np.nanmean(pi_along_loci)

            # pi_along_loci = list(map(lambda x: x.replace(0, np.NaN), pi_along_loci))

            # Compute average genome wide pi excluding 0 pi sites

            for index, item in enumerate(pi_along_loci):
                if item == 0.0:
                    pi_along_loci[index] = np.NaN

            genomewide_pi_variable_sites = np.nanmean(pi_along_loci)

            mean_depth = np.nanmean(depth_loci_sample1)
            n_total_loci = len(pi_along_loci)

            file_suffix_name = "Summary_Stats_Sidekick"

            summary_tests_pd = pd.DataFrame(
                np.nan,
                index=[0],
                columns=[
                    "Genomewide_pi",
                    "Genomewide_pi_variable_sites",
                    "Mean_depth",
                    "n_total_loci",
                ],
            )

            summary_tests_pd.iloc[:, 0] = genomewide_pi
            summary_tests_pd.iloc[:, 1] = genomewide_pi_variable_sites
            summary_tests_pd.iloc[:, 2] = mean_depth
            summary_tests_pd.iloc[:, 3] = n_total_loci

            summary_tests_pd.to_csv(
                output_dir
                + "Strain_"
                + strain
                + "_SampleID1_"
                + sample_id
                + "_Pi_"
                + calc_type
                + "_"
                + str(end_index)
                + "_"
                + file_suffix_name
                + ".csv.bz2",
                compression="bz2",
            )

        else:

            (
                loci_index,
                loci_name,
                gene_vec,
                pi_along_loci,
                depth_loci_sample1,
                depth_loci_sample2,
                freq_loci_sample1,
                freq_loci_sample2,
            ) = calculator_functions.compute_stat_for_sample_pair(
                depth_path=depth_path,
                freq_path=freq_path,
                snpinfo_path=snpinfo_path,
                gene_path=gene_path,
                pangenome_path=pangenome_path,
                blacklist_path=blacklist_path,
                sample_id_vec=input_sample_id_vec,
                end_index=end_index,
                calc_type=calc_type,
                statistic=calc_statistic,
                coverage_min=cov_min,
                summarymode=summarymode,
            )

            if pi_along_loci[0] == "Nothing":
                exit()
            print("writing to table")
            tests_pd = pd.DataFrame(
                np.nan,
                index=list(range(len(pi_along_loci))),
                columns=[
                    "loci_index",
                    "loci_name",
                    "gene_id",
                    "pi_along_loci",
                    "depth_loci_sample1",
                    "depth_loci_sample2",
                    "freq_loci_sample1",
                    "freq_loci_sample2",
                ],
            )

            tests_pd.iloc[:, 0] = loci_index
            tests_pd.iloc[:, 1] = loci_name
            tests_pd.iloc[:, 2] = gene_vec
            tests_pd.iloc[:, 3] = pi_along_loci
            tests_pd.iloc[:, 4] = depth_loci_sample1
            tests_pd.iloc[:, 5] = depth_loci_sample2
            tests_pd.iloc[:, 6] = freq_loci_sample1
            tests_pd.iloc[:, 7] = freq_loci_sample2

            genomewide_pi = np.nanmean(pi_along_loci)

            # pi_along_loci = list(map(lambda x: x.replace(0, np.NaN), pi_along_loci))

            # Compute average genome wide pi excluding 0 pi sites

            for index, item in enumerate(pi_along_loci):
                if item == 0.0:
                    pi_along_loci[index] = np.NaN

            genomewide_pi_variable_sites = np.nanmean(pi_along_loci)

            mean_depth1 = np.nanmean(depth_loci_sample1)
            mean_depth2 = np.nanmean(depth_loci_sample2)
            n_total_loci = len(pi_along_loci)

            file_suffix_name = "Summary_Stats_Sidekick"

            summary_tests_pd = pd.DataFrame(
                np.nan,
                index=[0],
                columns=[
                    "Genomewide_pi",
                    "Genomewide_pi_variable_sites",
                    "Mean_depth1",
                    "Mean_depth2",
                    "n_total_loci",
                ],
            )

            summary_tests_pd.iloc[:, 0] = genomewide_pi
            summary_tests_pd.iloc[:, 1] = genomewide_pi_variable_sites
            summary_tests_pd.iloc[:, 2] = mean_depth1
            summary_tests_pd.iloc[:, 3] = mean_depth2
            summary_tests_pd.iloc[:, 4] = n_total_loci

            summary_tests_pd.to_csv(
                output_dir
                + "Strain_"
                + strain
                + "_SampleID1_"
                + sample_id
                + "_SampleID2_"
                + sample_id2
                + "_Pi_"
                + calc_type
                + "_"
                + str(end_index)
                + "_"
                + file_suffix_name
                + ".csv.bz2",
                compression="bz2",
            )
        if loci_stats == True:
            if sample_id2 == "":

                tests_pd.to_csv(
                    output_dir
                    + "Strain_"
                    + strain
                    + "_SampleID1_"
                    + sample_id
                    + "_SinglePi_"
                    + calc_type
                    + "_"
                    + str(end_index)
                    + "_Loci_Stats.csv.bz2",
                    compression="bz2",
                )
            else:

                tests_pd.to_csv(
                    output_dir
                    + "Strain_"
                    + strain
                    + "_SampleID1_"
                    + sample_id
                    + "_SampleID2_"
                    + sample_id2
                    + "_PairedPi_"
                    + calc_type
                    + "_"
                    + str(end_index)
                    + "_Loci_Stats.csv.bz2",
                    compression="bz2",
                )

    elif calc_statistic == "Fst":

        if summarymode:
            (
                Fst_ratio_of_means,
                Fst_mean_of_ratios,
                n_total_loci,
            ) = calculator_functions.compute_stat_for_sample_pair(
                depth_path=depth_path,
                freq_path=freq_path,
                snpinfo_path=snpinfo_path,
                gene_path=gene_path,
                pangenome_path=pangenome_path,
                blacklist_path=blacklist_path,
                sample_id_vec=input_sample_id_vec,
                end_index=end_index,
                calc_type=calc_type,
                statistic=calc_statistic,
                coverage_min=cov_min,
                summarymode=summarymode,
            )
            if Fst_ratio_of_means == "Nothing":
                exit()

            print("mean of ratios")
            print(Fst_mean_of_ratios)

            print("ratio of means")
            print(Fst_ratio_of_means)

        elif calc_type in ["Regular", "Schloissnig"]:

            (
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
            ) = calculator_functions.compute_stat_for_sample_pair(
                depth_path=depth_path,
                freq_path=freq_path,
                snpinfo_path=snpinfo_path,
                gene_path=gene_path,
                pangenome_path=pangenome_path,
                blacklist_path=blacklist_path,
                sample_id_vec=input_sample_id_vec,
                end_index=end_index,
                calc_type=calc_type,
                statistic=calc_statistic,
                coverage_min=cov_min,
                summarymode=summarymode,
            )
            if fst_along_loci == "Nothing":
                exit()

            print("mean of fst_along_loci")
            print(np.nanmean(fst_along_loci))
            print("fst_across_loci from mean pi")
            print(fst_across_loci)

            tests_pd = pd.DataFrame(
                np.nan,
                index=list(range(len(fst_along_loci))),
                columns=[
                    "loci_index",
                    "loci_name",
                    "copy_number",
                    "gene_id",
                    "pi_within_along_loci",
                    "pi_between_along_loci",
                    "fst_along_loci",
                    "fst_across_loci",
                    "depth_loci_sample1",
                    "depth_loci_sample2",
                    "freq_loci_sample1",
                    "freq_loci_sample2",
                ],
            )

            start_data_col = 4
            tests_pd.iloc[:, 0] = loci_index
            tests_pd.iloc[:, 1] = loci_name

            tests_pd.iloc[:, 2] = cn_vec
            tests_pd.iloc[:, 3] = gene_vec

            pi_within_along_loci = [round(p_i, 5) for p_i in pi_within_along_loci]
            pi_between_along_loci = [round(p_i, 5) for p_i in pi_between_along_loci]
            fst_along_loci = [round(p_i, 5) for p_i in fst_along_loci]

            tests_pd.iloc[:, (0 + start_data_col)] = pi_within_along_loci
            tests_pd.iloc[:, (1 + start_data_col)] = pi_between_along_loci
            tests_pd.iloc[:, (2 + start_data_col)] = fst_along_loci
            tests_pd.iloc[:, (3 + start_data_col)] = list(
                np.repeat(round(fst_across_loci, 5), len(fst_along_loci))
            )
            tests_pd.iloc[:, (4 + start_data_col)] = depth_loci_sample1
            tests_pd.iloc[:, (5 + start_data_col)] = depth_loci_sample2
            tests_pd.iloc[:, (6 + start_data_col)] = freq_loci_sample1
            tests_pd.iloc[:, (7 + start_data_col)] = freq_loci_sample2

        elif calc_type == "Hudson":
            (
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
            ) = calculator_functions.compute_stat_for_sample_pair(
                depth_path=depth_path,
                freq_path=freq_path,
                snpinfo_path=snpinfo_path,
                gene_path=gene_path,
                pangenome_path=pangenome_path,
                blacklist_path=blacklist_path,
                sample_id_vec=input_sample_id_vec,
                end_index=end_index,
                calc_type=calc_type,
                statistic=calc_statistic,
                coverage_min=cov_min,
                summarymode=summarymode,
            )
            tests_pd = pd.DataFrame(
                np.nan,
                index=list(range(len(fst_along_loci))),
                columns=[
                    "loci_index",
                    "loci_name",
                    "copy_number",
                    "gene_id",
                    "fst_along_loci",
                    "fst_across_loci",
                    "depth_loci_sample1",
                    "depth_loci_sample2",
                    "freq_loci_sample1",
                    "freq_loci_sample2",
                    "pi_sample1",
                    "pi_sample2",
                ],
            )

            if fst_across_loci == "Nothing":
                exit()
            start_data_col = 4
            fst_along_loci = [round(p_i, 5) for p_i in fst_along_loci]
            tests_pd.iloc[:, 0] = loci_index
            tests_pd.iloc[:, 1] = loci_name
            tests_pd.iloc[:, 2] = cn_vec
            tests_pd.iloc[:, 3] = gene_vec

            tests_pd.iloc[:, (0 + start_data_col)] = fst_along_loci
            tests_pd.iloc[:, (1 + start_data_col)] = list(
                np.repeat(round(fst_across_loci, 5), len(fst_along_loci))
            )
            tests_pd.iloc[:, (2 + start_data_col)] = depth_loci_sample1
            tests_pd.iloc[:, (3 + start_data_col)] = depth_loci_sample2
            tests_pd.iloc[:, (4 + start_data_col)] = freq_loci_sample1
            tests_pd.iloc[:, (5 + start_data_col)] = freq_loci_sample2
            tests_pd.iloc[:, (6 + start_data_col)] = pi_sample1
            tests_pd.iloc[:, (7 + start_data_col)] = pi_sample2

        if summarymode:
            file_suffix_name = "Summary_Stats"

            summary_tests_pd = pd.DataFrame(
                np.nan,
                index=[0],
                columns=[" Fst_ratio_of_mean", "Fst_mean_of_ratios", "n_total_loci"],
            )

            summary_tests_pd.iloc[:, 0] = Fst_ratio_of_means
            summary_tests_pd.iloc[:, 1] = Fst_mean_of_ratios
            summary_tests_pd.iloc[:, 2] = n_total_loci

        else:
            file_suffix_name = "Loci_Stats"
            tests_pd.to_csv(
                output_dir
                + "Strain_"
                + strain
                + "_SampleID1_"
                + sample_id
                + "_SampleID2_"
                + sample_id2
                + "_Fst_"
                + calc_type
                + "_"
                + str(end_index)
                + "_"
                + file_suffix_name
                + ".csv.bz2",
                compression="bz2",
            )
            # summary stats
            file_suffix_name = "Summary_Stats_Sidekick"
            Fst_ratio_of_means = fst_across_loci
            Fst_mean_of_ratios = np.nanmean(fst_along_loci)
            n_total_loci = len(fst_along_loci)

            if calc_type == "Hudson":
                mean_pi_sample1 = np.nanmean(pi_sample1)
                mean_pi_sample2 = np.nanmean(pi_sample2)
                summary_tests_pd = pd.DataFrame(
                    np.nan,
                    index=[0],
                    columns=[
                        " Fst_ratio_of_mean",
                        "Fst_mean_of_ratios",
                        "n_total_loci",
                        "PiSample1",
                        "PiSample2",
                    ],
                )
                summary_tests_pd.iloc[:, 3] = mean_pi_sample1 #MW 08/28/24: this used to only be in the else portion of the ifelse statement
                summary_tests_pd.iloc[:, 4] = mean_pi_sample2 #MW 08/28/24: this used to only be in the else portion of the ifelse statement
            else:
                summary_tests_pd = pd.DataFrame(
                    np.nan,
                    index=[0],
                    columns=[
                        " Fst_ratio_of_mean",
                        "Fst_mean_of_ratios",
                        "n_total_loci",
                    ],
                )


            summary_tests_pd.iloc[:, 0] = Fst_ratio_of_means
            summary_tests_pd.iloc[:, 1] = Fst_mean_of_ratios
            summary_tests_pd.iloc[:, 2] = n_total_loci

        summary_tests_pd.to_csv(
            output_dir
            + "Strain_"
            + strain
            + "_SampleID1_"
            + sample_id
            + "_SampleID2_"
            + sample_id2
            + "_Fst_"
            + calc_type
            + "_"
            + str(end_index)
            + "_"
            + file_suffix_name
            + ".csv.bz2",
            compression="bz2",
        )

    print("done")
    print(output_dir)

    print(datetime.now() - start)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", help="sample id", type=str, default="")
    parser.add_argument(
        "--sample_id2", help="Minimal tester mode? (0 or 1)", type=str, default=""
    )
    parser.add_argument(
        "--end_index", help="Minimum read depth for a valid allele frequency", type=int
    )
    parser.add_argument("--strain", help="Strain to study", type=str)
    parser.add_argument("--study", help="Study to study", type=str)
    parser.add_argument("--calc_type", help="Method to calculate Pi and FST", type=str)
    parser.add_argument("--statistic", help="Pi or Fst", type=str)
    parser.add_argument("--cov_min", help="Coverage minimum", type=int)
    parser.add_argument(
        "--custom",
        help="Minimal tester mode? (0 or 1)",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--local",
        help="Minimal tester mode? (0 or 1)",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "--summarymode",
        help="Only include mean ratios and ratio means",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--loci_stats", 
        help="Compute loci-level statistics as well as summary.",
        default=False,
        action="store_true",
    )
                       

    args = parser.parse_args()
    sample_id = args.sample_id
    sample_id2 = args.sample_id2
    end_index = args.end_index
    study = args.study
    strain = args.strain
    calc_type = args.calc_type
    calc_statistic = args.statistic
    cov_min = args.cov_min
    custom = args.custom
    local = args.local
    summarymode = args.summarymode
    loci_stats = args.loci_stats

    main()

