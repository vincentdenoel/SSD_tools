import numpy as np
import pymc3 as pm

# Define the true parameters of the system
k_true = 1.0  # spring stiffness
c_true = 0.1  # damping coefficient

# Generate some simulated data
np.random.seed(42)  # for reproducibility
n_samples = 200
t = np.linspace(0, 10, n_samples)
x0 = 0.1
v0 = 0.0
w = np.random.normal(scale=0.1, size=n_samples)
x = x0*np.cos(np.sqrt(k_true)*t)*np.exp(-c_true*t/2) + v0*np.sin(np.sqrt(k_true)*t)*np.exp(-c_true*t/2) + w

# Define the prior distributions for the parameters
with pm.Model() as model:
    k = pm.Normal('k', mu=0.5, sd=0.5)
    c = pm.Normal('c', mu=0.05, sd=0.05)

    # Define the likelihood function
    likelihood = pm.Normal('likelihood', mu=x, sd=0.1, observed=x)

    # Run the MCMC algorithm to obtain posterior samples
    trace = pm.sample(1000, tune=1000, chains=2)

# Print the summary statistics of the posterior distributions
pm.summary(trace)
