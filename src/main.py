import os
import sys
import getopt

import parse_genet
import mcmc
import gigrnd


# Hardcoded parameters
ref_dir = "../data/1kg/ldblk_1kg_eur"
bim_prefix = "../data/test_data/test"
sst_file = "../data/test_data/sumstats_se.txt"
chrom = 22
n_gwas = 200000
out_dir = "../data/output_estimates/1kg/1kg_eur/"

# Default parameters
a = 1
b = 0.5
phi = None
n_iter = 1000
n_burnin = 500
thin = 5
beta_std = 'False'
write_psi = 'False'
seed = None


def main():
    print('##### process chromosome %d #####' % chrom)
    ref_dir_basename = os.path.basename(ref_dir.rstrip('/'))  # Strip the trailing slash if present

    if '1kg' in ref_dir_basename.lower():  # Using .lower() to make it case-insensitive
       ref_dict = parse_genet.parse_ref(ref_dir + '/snpinfo_1kg_hm3', chrom)
    elif 'ukbb' in ref_dir_basename.lower():  # Using .lower() to make it case-insensitive
       ref_dict = parse_genet.parse_ref(ref_dir + '/snpinfo_ukbb_hm3', chrom)
    else:
       raise ValueError(f"Unable to determine reference dictionary for chromosome {chrom} with ref_dir {ref_dir}")

    if ref_dict is None:
       raise ValueError(f"Reference dictionary not created for chromosome {chrom}. Check the reference directory name and format.")

    vld_dict = parse_genet.parse_bim(bim_prefix, chrom)
    sst_dict = parse_genet.parse_sumstats(ref_dict, vld_dict, sst_file, n_gwas)
    ld_blk, blk_size = parse_genet.parse_ldblk(ref_dir, sst_dict, chrom)

    mcmc.mcmc(a, b, phi, sst_dict, n_gwas, ld_blk, blk_size, n_iter, n_burnin, thin, chrom, out_dir, beta_std, write_psi, seed)
    print('\n')


if __name__ == '__main__':
    main()
