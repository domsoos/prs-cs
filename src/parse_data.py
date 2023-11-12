import os
import scipy as sp
from scipy.stats import norm
from scipy import linalg
import h5py

class GeneticDataParser:
    def __init__(self, chromosome):
        self.chromosome = chromosome
        self.ATGC = ['A', 'T', 'G', 'C']
        self.mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def parse_reference_file(self, ref_file):
        ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[]}
        with open(ref_file) as file:
            next(file)  # Skipping header
            for line in file:
                ll = line.strip().split()
                if int(ll[0]) == self.chromosome:
                    self._add_ref_data(ref_dict, ll)
        return ref_dict

    def _add_ref_data(self, ref_dict, ll):
        ref_dict['CHR'].append(self.chromosome)
        ref_dict['SNP'].append(ll[1])
        ref_dict['BP'].append(int(ll[2]))
        ref_dict['A1'].append(ll[3])
        ref_dict['A2'].append(ll[4])
        ref_dict['MAF'].append(float(ll[5]))

    def parse_bim(self, bim_file):
        vld_dict = {'SNP':[], 'A1':[], 'A2':[]}
        with open(bim_file + '.bim') as file:
            for line in file:
                ll = line.strip().split()
                if int(ll[0]) == self.chromosome:
                    vld_dict['SNP'].append(ll[1])
                    vld_dict['A1'].append(ll[4])
                    vld_dict['A2'].append(ll[5])
        return vld_dict

    def parse_sumstats(self, ref_dict, vld_dict, sst_file, n_subj):
        sst_dict = {'SNP':[], 'A1':[], 'A2':[]}
        with open(sst_file) as file:
            next(file)  # Skipping header
            for line in file:
                ll = line.strip().split()
                if ll[1] in self.ATGC and ll[2] in self.ATGC:
                    sst_dict['SNP'].append(ll[0])
                    sst_dict['A1'].append(ll[1])
                    sst_dict['A2'].append(ll[2])

        vld_snp, ref_snp, sst_snp, comm_snp = self._get_snp_sets(vld_dict, ref_dict, sst_dict)

        sst_eff, sst_dict_final = self._compute_eff(sst_file, comm_snp, n_subj, header)
        for (ii, snp) in enumerate(ref_dict['SNP']):
            if snp in sst_eff:
                self._add_sst_data(sst_dict_final, ref_dict, sst_eff, ii, snp, comm_snp)
        return sst_dict_final