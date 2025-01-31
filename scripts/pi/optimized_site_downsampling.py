import argparse
import os
import bz2
import sys
import random
import numpy as np
from datetime import datetime
import shutil
def select_mice(mice_numbers, colnames):
    mice_ID = ["M" + str(m) for m in mice_numbers]
    return [i for i, col in enumerate(colnames) if any(m in col for m in mice_ID)]

def generate_reads_vector(target_depth, target_freq):
    number_ref = int(round(target_depth * target_freq))
    reads_vector = np.concatenate([np.ones(number_ref, dtype=int), np.zeros(target_depth - number_ref, dtype=int)])
    np.random.shuffle(reads_vector)
    return reads_vector

def generate_new_depth_and_freq(simulated_reads, median_depth):
    sampled_reads = np.random.choice(simulated_reads, median_depth, replace=False)
    out_freq = np.mean(sampled_reads)
    return len(sampled_reads), out_freq

def main():
    depth_path = os.path.join(in_data_path, strain, 'snps_depth.txt.bz2')
    freq_path = os.path.join(in_data_path, strain, 'snps_ref_freq.txt.bz2')

    depth_out_path = os.path.join(out_data_path, strain, 'snps_depth.txt.bz2')
    freq_out_path = os.path.join(out_data_path, strain, 'snps_ref_freq.txt.bz2')

    os.makedirs(os.path.join(out_data_path, strain), exist_ok=True)

    with bz2.open(depth_path, 'rt') as depth:
        headers = depth.readline().strip().split('\t')

    select_indices = select_mice([1, 2, 3, 4, 5, 6], headers)
    select_indices = [s-1 for s in select_indices]
    target_indices = select_mice([7, 8], headers)
    target_indices = [s-1 for s in target_indices]

    if len(select_indices)==0 or len(target_indices)  == 0:
        shutil.copy(depth_path, depth_out_path)
        shutil.copy(freq_path, freq_out_path)
        sys.exit()




    measured_sites = 0
    total_sites = 0
    start = datetime.now()

    with bz2.open(depth_path, 'rt') as depth, bz2.open(freq_path, 'rt') as freq, \
         bz2.open(depth_out_path, 'wt') as depth_out, bz2.open(freq_out_path, 'wt') as freq_out:

        depth_out.write('\t'.join(headers) + '\n')
        freq_out.write('\t'.join(headers) + '\n')
        next(depth)  # Skip header
        next(freq)

        for line_depth, line_freq in zip(depth, freq):
            total_sites += 1
            if total_sites % 10000 == 0:
                print(f"Elapsed time: {datetime.now() - start}, Measured sites: {measured_sites}, Total sites: {total_sites}")

            # Process depth line
            line_clean = line_depth.strip().split('\t')
            site_name = line_clean[0]
            values_depth = np.array(list(map(int, line_clean[1:])), dtype=int)
            select_depth = values_depth[select_indices]
            if metric == "median":
                median_depth = int(np.median(select_depth))
            else:
                median_depth = round(np.mean(select_depth))
                # print(np.mean(select_depth))
                # print(int(np.mean(select_depth)))

            zero_out_flag = median_depth < 4
            if not zero_out_flag:
                measured_sites += 1

            # Process frequency line
            values_freq = np.array(list(map(float, line_freq.strip().split('\t')[1:])), dtype=float)
            new_values_freq = values_freq.copy()
            new_values_depth = values_depth.copy()

            if zero_out_flag:
                new_values_freq[target_indices] = 0
                new_values_depth[target_indices] = 0
            else:
                for idx in target_indices:
                    target_freq = values_freq[idx]
                    target_depth = values_depth[idx]
                    if target_depth > median_depth:
                        simulated_reads = generate_reads_vector(int(target_depth), target_freq)
                        new_depth, new_freq = generate_new_depth_and_freq(simulated_reads, median_depth)
                        new_values_freq[idx] = new_freq
                        new_values_depth[idx] = new_depth
                #check = 0

            depth_out.write(f"{site_name}\t" + '\t'.join(map(str, new_values_depth)) + '\n')
            freq_out.write(f"{site_name}\t" + '\t'.join(map(str, new_values_freq)) + '\n')

    print(f"Measured sites: {measured_sites}")
    print(f"Original sites: {total_sites}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--strain", help="strain", type=str, default="")
    parser.add_argument("--metric", help="median", type=str, default="")
    parser.add_argument("--local", help="strain", action='store_true', default=False)
    
    
    args = parser.parse_args()
    strain = args.strain
    metric = args.metric
    local = args.local
   
    in_data_path = ~/merged_data/snps/
    if metric == "median":
        out_data_path = "~/SiteDownsampled/merged_data/snps"
    else: 
        out_data_path = "~/SiteDownsampledByMean/merged_data/snps"

  
    if not os.path.exists(out_data_path):
        os.makedirs(out_data_path)
    main()
