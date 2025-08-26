import dadi
import numpy as np

# 1. Load your folded SFS
# load your counts vector
counts = np.loadtxt('BP1.folded.sfs')   # shape (31,)

# construct a 1D Spectrum
data = dadi.Spectrum(counts)

print("Loaded SFS:", data)
print("Sample size (ns):", data.sample_sizes)

ns=data.sample_sizes  # sample sizes tuple

# 2. Define a simple bottleneck model
def bottleneck_model(params, ns, pts):
    """
    params = (nuB, TB, nuF, TF)
      nuB = Ne during bottleneck (relative to ancestral Ne)
      TB  = duration of bottleneck (in 2*Na generations)
      nuF = final size after recovery (relative to ancestral Ne)
      TF  = duration of recovery (in 2*Na generations)
    ns   = sample sizes tuple
    pts  = grid point count for extrapolation
    """
    nuB, TB, nuF, TF = params
    # grid for extrapolation
    xx = dadi.Numerics.default_grid(pts)
    # ancestral equilibrium
    phi = dadi.PhiManip.phi_1D(xx)            # <— use phi_1D for equilibrium
    # integrate through bottleneck
    phi = dadi.Integration.one_pop(phi, xx, TB, nu=nuB)
    # integrate through recovery
    phi = dadi.Integration.one_pop(phi, xx, TF, nu=nuF)
    # get model SFS
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

# 3. Set up grid & initial guesses
pts_l = [max(ns)*2 + 10]        # e.g. [70] for smooth extrapolation
p0   = [0.1, 0.05, 1.0, 0.1]     # initial: 10% Ne during bottleneck, bottleneck 0.05, recover to 1×, recover time 0.1
lower_bounds = [1e-3, 1e-4, 0.5, 1e-4]
upper_bounds = [10,    1,    10,   1]

# 4. Optimize parameters to fit data
func_ex = dadi.Numerics.make_extrap_log_func(bottleneck_model)
popt = dadi.Inference.optimize_log(
    p0, data, func_ex, pts_l,
    lower_bound=lower_bounds,
    upper_bound=upper_bounds,
    verbose=1, maxiter=100
)
print(f"\nOptimized parameters: {popt}")

# 5. Compare model vs. data
model = func_ex(popt, ns, pts_l)
ll_model = dadi.Inference.ll_multinom(model, data)
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print(f"Log‐likelihood = {ll_model:.2f}, θ = {theta:.2f}")

# 6. Save a simple plot (optional)
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6,4))
# basic comparison plot
dadi.Plotting.plot_1d_comp_multinom(model, data)

# Optional: add a title or axis labels
plt.title("Bottleneck fit for BP1")
plt.xlabel("Minor allele count")
plt.ylabel("Frequency")

plt.tight_layout()
plt.savefig('bottleneck_fit_BP1.png', dpi=300)


