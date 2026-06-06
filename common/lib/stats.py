"""MD trajectory statistics with autocorrelation correction.

Best-practice §32: MD trajectories are highly autocorrelated time series.
Do NOT use independent-sample t-tests. Report effective sample size (n_eff)
and use correlated t-tests or block-bootstrap for significance testing.
"""
import numpy as np
from typing import Sequence


def autocorr_time(x: np.ndarray, max_lag: int = None) -> float:
    """Estimate integrated autocorrelation time (IACT) using Geyer's method.

    Returns tau in units of *frames*.  tau = 1 means uncorrelated.
    """
    x = np.asarray(x)
    if x.size < 10:
        return 1.0
    x = x - x.mean()
    c0 = np.mean(x * x)
    if c0 == 0:
        return 1.0

    max_lag = max_lag or min(len(x) // 4, 1000)
    c = np.correlate(x, x, mode='full')[len(x) - 1:]
    c = c / c[0]  # normalise

    tau = 1.0
    for lag in range(1, max_lag, 2):
        if lag + 1 >= len(c):
            break
        gamma = c[lag] + c[lag + 1]
        if gamma <= 0:
            break
        tau += 2 * gamma
    return tau


def effective_sample_size(x: np.ndarray) -> float:
    """n_eff = n / (1 + 2 * sum_autocorr) ≈ n / tau"""
    x = np.asarray(x)
    if x.size < 2:
        return float(x.size)
    tau = autocorr_time(x)
    return max(1.0, len(x) / tau)


def summarize(x: np.ndarray, name: str = "") -> dict:
    """Return a dict with mean, std, median, IQR, n_eff for a 1-D time series.

    Use this instead of raw mean±std when reporting MD observables.
    """
    x = np.asarray(x)
    n = x.size
    if n == 0:
        return {"name": name, "n": 0, "mean": np.nan, "std": np.nan,
                "median": np.nan, "q25": np.nan, "q75": np.nan,
                "iqr": np.nan, "tau_frames": np.nan, "n_eff": np.nan}

    tau = autocorr_time(x)
    n_eff = max(1.0, n / tau)
    return {
        "name": name,
        "n": n,
        "mean": float(x.mean()),
        "std": float(x.std(ddof=1)),
        "median": float(np.median(x)),
        "q25": float(np.percentile(x, 25)),
        "q75": float(np.percentile(x, 75)),
        "iqr": float(np.percentile(x, 75) - np.percentile(x, 25)),
        "tau_frames": float(tau),
        "n_eff": float(n_eff),
    }


def replica_cv(replica_means: Sequence[float]) -> float:
    """Coefficient of variation across replica means.

    Best-practice §11: report replica CV, not just within-replica std.
    """
    arr = np.asarray(replica_means, dtype=float)
    if len(arr) < 2 or arr.mean() == 0:
        return np.nan
    return float(arr.std(ddof=1) / abs(arr.mean()))


def format_summary(s: dict, precision: int = 2) -> str:
    """Pretty-print a summary dict."""
    p = precision
    return (f"{s['name']:<20s}  "
            f"mean={s['mean']:.{p}f}  std={s['std']:.{p}f}  "
            f"median={s['median']:.{p}f}  IQR=[{s['q25']:.{p}f},{s['q75']:.{p}f}]  "
            f"n_eff={s['n_eff']:.0f}")


def correlated_t_test(x: np.ndarray, y: np.ndarray) -> dict:
    """Paired correlated t-test for two MD trajectories of the same length.

    Delong et al. (2011) / best-practice §32.
    Returns t-statistic, df, and p-value (two-tailed).
    """
    x, y = np.asarray(x), np.asarray(y)
    if len(x) != len(y):
        raise ValueError("x and y must have same length")
    d = x - y
    n = len(d)
    if n < 2:
        return {"t": np.nan, "df": np.nan, "p": np.nan}

    d_mean = d.mean()
    tau = autocorr_time(d)
    n_eff = max(1.0, n / tau)
    # Use effective sample size for std error
    se = d.std(ddof=1) / np.sqrt(n_eff)
    if se == 0:
        return {"t": np.nan, "df": np.nan, "p": np.nan}

    t_stat = d_mean / se
    df = n_eff - 1
    from scipy import stats as st
    p_val = 2 * st.t.cdf(-abs(t_stat), df)
    return {"t": float(t_stat), "df": float(df), "p": float(p_val), "n_eff": float(n_eff)}
