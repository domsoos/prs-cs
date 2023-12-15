import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t, norm, cauchy
from scipy.special import gamma

# Define the range for tau and beta
tau_range = np.linspace(0, 1, 1000)
beta_range = np.linspace(-3, 3, 1000)
beta_tail_range = np.linspace(3, 7, 1000)

# Define the parameters for the priors
a_values = [1/2, 1, 3/2]  # Values of a in the prior density functions
phi = 1  # Global shrinkage parameter phi is set to 1
b = 1/2  # Value of b in the prior density function

# Function to calculate the Three-Parameter Beta prior density
def tpb_density(tau, a, b, phi):
    return (gamma(a + b) / (gamma(a) * gamma(b))) * \
           (phi**b) * (tau**(b - 1)) * ((1 - tau)**(a - 1)) * \
           (1 + (phi - 1) * tau)**(-(a + b))

# Define the prior density function for beta in the central region
def beta_density_central(beta, a, phi):
    # Using the Cauchy distribution as an example for the horseshoe prior
    if a == 1/2:
        return cauchy.pdf(beta)
    else:
        return t.pdf(beta, df=2 * a, scale=np.sqrt(1 / (a * phi)))
def prior_density_beta(beta, a, b, phi):
    # Since a and b are scalars, we can use the scalar version of the cauchy.pdf
    return cauchy.pdf(beta) * tpb_density((1/(1 + phi))**a, a, b, phi)

# Define the prior density function for beta in the tails
def beta_density_tails(beta, a, phi):
    return t.pdf(beta, df=2 * a, scale=np.sqrt(1 / (a * phi)))

# Initialize plots
fig, axs = plt.subplots(3, 1, figsize=(6, 10))

# Plot for prior density of τ_j (shrinkage factor)
for a in a_values:
    tau_density = tpb_density(tau_range, a, b, phi)
    axs[0].plot(tau_range, tau_density, label=f'a = {a}')

axs[0].set_title('Prior density of $\\tau_j$')
axs[0].legend()

# Plot for prior density of β_j in the central region
for a in a_values:
    beta_central_density = beta_density_central(beta_range, a, phi)
    axs[1].plot(beta_range, beta_central_density, label=f'a = {a}')

# Adding standard normal for comparison
axs[1].plot(beta_range, norm.pdf(beta_range), 'k--', label='Normal')
axs[1].set_title('Prior density of $\\beta_j$: central region')
axs[1].legend()

# Plot for prior density of β_j in the tails
for a in a_values[::-1]:  # Reversing the order for correct layering
    beta_tail_density = beta_density_tails(beta_tail_range, a, phi)
    axs[2].plot(beta_tail_range, beta_tail_density, label=f'a = {a}')

axs[2].set_title('Prior density of $\\beta_j$: tails')
axs[2].legend()

plt.tight_layout()
plt.show()
