import numpy as np
import pymc as pm
import matplotlib.pyplot as plt

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

# Define the model
with pm.Model() as model:
    # Priors for the parameters
    k = pm.Normal('k', mu=1.0, sigma=0.5)
    c = pm.Normal('c', mu=0.1, sigma=0.05)

    # Expected value of outcome
    x_pred = x0 * pm.math.cos(pm.math.sqrt(k) * t) * pm.math.exp(-c * t / 2) + v0 * pm.math.sin(pm.math.sqrt(k) * t) * pm.math.exp(-c * t / 2)

    # Likelihood (sampling distribution) of observations
    likelihood = pm.Normal('y', mu=x_pred, sigma=0.1, observed=x)

# Print the summary statistics of the posterior distributions
with model:
    trace = pm.sample(2000, tune=1000, return_inferencedata=False)
    pm.summary(trace).round(2)
# Plot the posterior distributions of the parameters
pm.plot_posterior(trace, var_names=['k', 'c'], hdi_prob=0.95)
# Show the plot

plt.show() 

print('Done.')
