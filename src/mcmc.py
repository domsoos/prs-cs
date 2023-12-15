import scipy as sp
from scipy import linalg
from numpy import random
import numpy as np
import gigrnd

def update_shrinkage_params(delta_param, shrinkage_param, num_snps, prior_shape_param, sample_size, effect_size, variance_param):
    """
    Update shrinkage parameters based on the gamma-gamma prior.
    :param delta_param: Gamma distribution parameters for updating shrinkage params.
    :param shrinkage_param_estimated: Current estimate of the shrinkage parameters.
    :param num_snps: Number of SNPs.
    :param prior_shape_param: Shape parameter for the prior.
    :param sample_size: Number of samples.
    :param effect_size_estimated: Current estimate of the effect sizes.
    :param variance_param_estimated: Current estimate of the variance.
    :return: Updated shrinkage parameters.
    """
    for jj in range(num_snps):
        shrinkage_param[jj] = gigrnd.gigrnd(prior_shape_param - 0.5, 2.0 * delta_param[jj], sample_size * effect_size[jj] ** 2 / variance_param)
    shrinkage_param[shrinkage_param>1] = 1.0
    return shrinkage_param

def update_effect_sizes(effect_size, ld_blocks, block_size, num_blocks, shrinkage_param, effect_size_marginal, variance_param, sample_size):
    """
    Update effect sizes using the LD matrix and shrinkage parameters.
    :param ld_blocks: LD matrix blocks.
    :param block_size: Sizes of LD blocks.
    :param num_blocks: Number of LD blocks.
    :param shrinkage_param_estimated: Current shrinkage parameters.
    :param effect_size_marginal: Marginal effect sizes.
    :param variance_param_estimated: Current variance estimate.
    :param sample_size: Number of samples.
    :return: Updated effect sizes and the quadratic term for variance update.
    """
    quadratic_term = 0.0
    mm = 0
    for block_index in range(num_blocks):
        if block_size[block_index] == 0:
            continue
        block_range = range(mm, mm + block_size[block_index])
        ld_block_inv = ld_blocks[block_index] + sp.diag(1.0 / shrinkage_param[block_range].T[0])
        cholesky_factor = linalg.cholesky(ld_block_inv)
        temporary_effect_size = linalg.solve_triangular(cholesky_factor, effect_size_marginal[block_range], trans='T') + sp.sqrt(variance_param / sample_size) * random.randn(len(block_range), 1)
        effect_size[block_range] = linalg.solve_triangular(cholesky_factor, temporary_effect_size, trans='N')
        quadratic_term += sp.dot(sp.dot(effect_size[block_range].T, ld_block_inv), effect_size[block_range])
        mm += block_size[block_index]
    return effect_size, quadratic_term

# Function to calculate the shrinkage factor based on Equation 8
def calculate_shrinkage_factor(global_scaling_param, local_shrinkage_params):
    """
    Calculate shrinkage factor for each marker.
    :param global_scaling_param: Global scaling parameter (phi).
    :param local_shrinkage_params: Local shrinkage parameters (psi).
    :return: Shrinkage factor for each marker.
    """
    return 1 / (1 + global_scaling_param * local_shrinkage_params)

# Function for Gamma-Gamma Prior based on Equation 9
def gamma_gamma_prior(a, b, num_snps):
    """
    Calculate gamma-gamma prior for local shrinkage parameters.
    :param a: Shape parameter alpha.
    :param b: Scale parameter beta.
    :param num_snps: Number of SNPs.
    :return: Generated local shrinkage parameters.
    """
    scaled_prior = prior_shape_param + prior_scale_param
    dominator = (shrinkage_param + global_shrinkage_param)
    return np.random.gamma(scaled_prior, 1.0 / dominator)

# Function to calculate the posterior estimates
def calculate_posterior_estimates(effect_size, effect_size_estimated, shrinkage_param, shrinkage_param_estimated, variance_param, variance_param_estimated, global_shrinkage_param, global_shrinkage_param_estimated):
    effect_size_estimated = effect_size_estimated + effect_size / num_posterior_stats
    shrinkage_param_estimated = shrinkage_param_estimated + shrinkage_param / num_posterior_stats
    variance_param_estimated = variance_param_estimated + variance_param / num_posterior_stats
    global_shrinkage_param_estimated = global_shrinkage_param_estimated + global_shrinkage_param / num_posterior_stats
    return effect_size_estimated, shrinkage_param_estimated, variance_param_estimated, global_shrinkage_param_estimated

def mcmc(prior_shape_param, prior_scale_param, global_shrinkage_param, sum_stat,sample_size, ld_blocks, block_size, num_iterations, burn_in, thinning_factor, chromosome_num, output_dir, effect_size_standardized, save_shrinkage_params, random_seed):
    print('Starting Polygenic MCMC...')

    if random_seed is not None:
        random.seed(random_seed)

    effect_size_marginal = sp.array(sum_stat['BETA'], ndmin=2).T
    minor_allele_frequency = sp.array(sum_stat['MAF'], ndmin=2).T
    num_posterior_stats = (num_iterations - burn_in) // thinning_factor
    num_snps = len(sum_stat['SNP'])
    num_blocks = len(ld_blocks)

    effect_size = sp.zeros((num_snps,1))
    shrinkage_param = sp.ones((num_snps,1))
    variance_param = 1.0
    variance_param_estimated = 0.0
    if global_shrinkage_param is None:
        global_shrinkage_param = 1.0
        phi_updt = True
    else:
        phi_updt = False
    
    global_shrinkage_param_estimated = 0.0

    effect_size_estimated = sp.zeros((num_snps, 1))
    shrinkage_param_estimated = sp.zeros((num_snps, 1))

    for iteration in range(1, num_iterations + 1):
        if iteration % 100 == 0:
        	if iteration > 500:
        		print(f'--- iter-{str(iteration)} ---\nlocal shrinkage := {shrinkage_param_estimated[:3]}\nglobal shrinkage := {global_shrinkage_param_estimated}')
        	else:
        		print(f'--- iter-{str(iteration)}')
        		
        # Update effect sizes and get quadratic term
        effect_size_temp, quadratic_term = update_effect_sizes(effect_size, ld_blocks, block_size, num_blocks, shrinkage_param, effect_size_marginal, variance_param, sample_size)

        # Update the gamma-gamma parameters
        delta_param = np.random.gamma(prior_shape_param + prior_scale_param, 1.0 / (shrinkage_param + global_shrinkage_param))

        # Update shrinkage parameters
        shrinkage_param = update_shrinkage_params(delta_param, shrinkage_param, num_snps, prior_shape_param, sample_size, effect_size, variance_param)

        # Update variance parameter
        error_term = max(sample_size / 2.0 * (1.0 - 2.0 * sum(effect_size * effect_size_marginal) + quadratic_term), sample_size / 2.0 * sum(effect_size ** 2 / shrinkage_param))
        variance_param = 1.0 / random.gamma((sample_size + num_snps) / 2.0, 1.0 / error_term)

        if phi_updt:
            w_param = random.gamma(1.0, 1.0 / (global_shrinkage_param + 1.0))
            global_shrinkage_param = random.gamma(num_snps * prior_scale_param + 0.5, 1.0 / (sum(delta_param) + w_param))

        if (iteration > burn_in) and (iteration % thinning_factor == 0):
            effect_size_estimated = effect_size_estimated + effect_size / num_posterior_stats
            shrinkage_param_estimated = shrinkage_param_estimated + shrinkage_param / num_posterior_stats
            variance_param_estimated = variance_param_estimated + variance_param / num_posterior_stats
            global_shrinkage_param_estimated = global_shrinkage_param_estimated + global_shrinkage_param / num_posterior_stats

    if effect_size_standardized == 'False':
        effect_size_estimated /= sp.sqrt(2.0 * minor_allele_frequency * (1.0 - minor_allele_frequency))

    if phi_updt:
        print(f'Estimated global shrinkage parameter: {global_shrinkage_param_estimated}')

    print('MCMC Complete.')

