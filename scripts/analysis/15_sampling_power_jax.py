"""How much sampling does it take to resolve the MM1/MM2 orientation difference?

Usage: python 15_sampling_power_jax.py <analysis_dir> <out_dir>

Runs on the GPU through JAX. Three questions, in order:

1. How correlated is the tilt signal in time? The integrated autocorrelation time tau_int is found
   per replica by FFT autocorrelation with Sokal's automatic window (window at the first M where
   M >= c * tau, c = 6). The effective sample size of a trajectory of length T is then
   N_eff = T / (2 tau_int), not the number of frames written.

2. Do the replicas agree to within what that predicts? The within-replica standard error is
   sigma * sqrt(2 tau_int / T). If the observed scatter between replica means is much larger, a slow
   coordinate is not being sampled inside one replica, and adding frames to a single run cannot fix
   it. The ratio of the two is the number that matters.

3. Given both, how long must a run be to resolve a difference of a given size? Reported for a
   two-motor comparison with n replicas of length T each, using the variance that actually applies,
   requiring the separation to exceed 2 standard errors of the difference.

Everything is computed twice, once from within-replica statistics and once from the between-replica
scatter, because the gap between them is the finding.
"""
import json
import os
import sys

import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)

A, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
DT_NS = 0.01          # 10 ps between rows of series.dat
BURN_NS = 30.0        # analysis window used throughout the project
TARGETS = (45.0, 20.0, 8.8)   # degrees: the published contrast, a modest one, and what we measured


@jax.jit
def autocorr(x):
    """Normalised autocorrelation of a 1-D series, via FFT (Wiener-Khinchin)."""
    x = x - x.mean()
    n = x.shape[0]
    f = jnp.fft.rfft(x, n=2 * n)
    ac = jnp.fft.irfft(f * jnp.conj(f), n=2 * n)[:n]
    return ac / ac[0]


def tau_int(x, c=6.0):
    """Integrated autocorrelation time in frames, Sokal automatic windowing."""
    rho = np.asarray(autocorr(jnp.asarray(x)))
    csum = 2.0 * np.cumsum(rho) - 1.0          # tau(M) = 1 + 2 sum_{t=1..M} rho(t)
    m = np.arange(len(csum))
    ok = m >= c * np.maximum(csum, 1e-12)
    w = int(np.argmax(ok)) if ok.any() else len(csum) - 1
    return max(float(csum[w]), 1.0), w


rows, series = [], {}
for M in ("MM1", "MM2"):
    for r in range(1, 6):
        s = np.loadtxt(os.path.join(A, "series", f"{M}_r{r}_series.dat"))
        keep = s[:, 0] >= BURN_NS * 1000.0
        tilt = s[keep, 1]
        t_frames, win = tau_int(tilt)
        T_ns = len(tilt) * DT_NS
        tau_ns = t_frames * DT_NS
        rows.append({"motor": M, "replica": r, "T_ns": T_ns, "mean_deg": float(tilt.mean()),
                     "sd_deg": float(tilt.std(ddof=1)), "tau_int_ns": tau_ns,
                     "n_eff": T_ns / (2 * tau_ns), "window_frames": win,
                     "sem_within_deg": float(tilt.std(ddof=1)) * np.sqrt(2 * tau_ns / T_ns)})
        series.setdefault(M, []).append(tilt)

out = {"window_ns": [BURN_NS, 100.0], "dt_ns": DT_NS, "per_replica": rows, "per_motor": {},
       "device": str(jax.devices()[0])}

for M in ("MM1", "MM2"):
    rs = [x for x in rows if x["motor"] == M]
    means = np.array([x["mean_deg"] for x in rs])
    sd_within = float(np.mean([x["sd_deg"] for x in rs]))
    tau = float(np.mean([x["tau_int_ns"] for x in rs]))
    T = float(np.mean([x["T_ns"] for x in rs]))
    sem_within = sd_within * np.sqrt(2 * tau / T)          # predicted scatter of one replica mean
    sd_between = float(means.std(ddof=1))                  # observed scatter of replica means
    out["per_motor"][M] = {
        "replica_means_deg": means.round(2).tolist(),
        "mean_deg": float(means.mean()), "sem_over_replicas_deg": float(sd_between / np.sqrt(len(means))),
        "sd_within_replica_deg": sd_within, "tau_int_ns": tau, "n_eff_per_replica": T / (2 * tau),
        "sem_predicted_from_one_replica_deg": sem_within,
        "sd_between_replica_means_deg": sd_between,
        "underdispersion_ratio": sd_between / sem_within,
    }

# --- how long would it take? -------------------------------------------------------------------
# Two motors, n replicas each of length T. SE of the difference of the two grand means:
#   from within-replica statistics : sqrt(2) * sigma * sqrt(2 tau / (n T))
#   from between-replica scatter   : sqrt(2) * sd_between / sqrt(n)        (T-independent)
sd_w = float(np.mean([out["per_motor"][m]["sd_within_replica_deg"] for m in ("MM1", "MM2")]))
tau_m = float(np.mean([out["per_motor"][m]["tau_int_ns"] for m in ("MM1", "MM2")]))
sd_b = float(np.mean([out["per_motor"][m]["sd_between_replica_means_deg"] for m in ("MM1", "MM2")]))
req = []
for delta in TARGETS:
    for n in (5, 10, 20):
        # need 2 * SE(difference) <= delta
        T_need = 8.0 * tau_m * sd_w ** 2 / (n * (delta / 2.0) ** 2)
        n_need_between = 8.0 * sd_b ** 2 / delta ** 2
        req.append({"delta_deg": delta, "n_replicas": n,
                    "T_per_replica_ns_if_only_fast_motion": T_need,
                    "total_ns_if_only_fast_motion": n * T_need,
                    "n_replicas_needed_given_observed_scatter": n_need_between})
out["requirements"] = req
out["inputs_for_requirements"] = {"sd_within_deg": sd_w, "tau_int_ns": tau_m, "sd_between_deg": sd_b}

json.dump(out, open(os.path.join(OUT, "sampling-power.json"), "w"), indent=1)

print("device:", out["device"])
print("\nper replica (30-100 ns window)")
print("  motor rep   mean    sd   tau_int   N_eff   SEM(within)")
for x in rows:
    print("  %-5s %2d  %6.1f %5.1f  %6.2f ns %6.1f      %5.2f deg"
          % (x["motor"], x["replica"], x["mean_deg"], x["sd_deg"], x["tau_int_ns"], x["n_eff"], x["sem_within_deg"]))
print("\nper motor")
for M in ("MM1", "MM2"):
    d = out["per_motor"][M]
    print("  %s: tau_int %.2f ns -> N_eff %.0f per 70 ns replica" % (M, d["tau_int_ns"], d["n_eff_per_replica"]))
    print("      one replica's mean should scatter by %.2f deg; observed scatter is %.2f deg (%.1fx)"
          % (d["sem_predicted_from_one_replica_deg"], d["sd_between_replica_means_deg"], d["underdispersion_ratio"]))
print("\nto resolve a difference at 2 standard errors")
for x in req:
    if x["n_replicas"] == 5:
        print("  %.1f deg: %.0f ns per replica x 5 if only fast motion mattered (%.1f us total);"
              % (x["delta_deg"], x["T_per_replica_ns_if_only_fast_motion"], x["total_ns_if_only_fast_motion"] / 1000.0))
        print("            but given the observed replica scatter, %.0f independent replicas are needed"
              % x["n_replicas_needed_given_observed_scatter"])
