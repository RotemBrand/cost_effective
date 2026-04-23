import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB
from sklearn.cluster import kmeans_plusplus


def balanced_kmeans_gurobi(
        X,
        k,
        max_iter=20,
        tol=1e-4,
        random_state=0,
        chain_len_sigma=0.0,
    ):
    """
    Performs balanced k-means clustering using Gurobi to solve the assignment problem at each iteration.

    This algorithm clusters points into exactly k clusters, enforcing balanced cluster sizes (or controlled diversity via chain_len_sigma).
    At each iteration, it solves a transportation linear program to optimally assign points to clusters given current centers,
    then updates cluster centers as the mean of assigned points.

    Parameters
    ----------
    X : np.ndarray
        Array of shape (n_samples, n_features) representing the data points.
    k : int
        Number of clusters.
    max_iter : int, optional
        Maximum number of k-means iterations (default: 20).
    tol : float, optional
        Tolerance for convergence based on center shift (default: 1e-4).
    random_state : int, optional
        Random seed for reproducibility (default: 0).
    chain_len_sigma : float, optional
        Standard deviation for cluster size diversity (default: 0.0, perfectly balanced).

    Returns
    -------
    chains : np.ndarray
        Array of shape (n_samples,) with cluster labels for each point.
    centers : np.ndarray
        Array of shape (k, n_features) with final cluster centers.
    """
    X = np.asarray(X, dtype=np.float32)
    # normelize
    _min = X.min(axis=0)
    _X = X - _min
    _max = _X.max(axis=0)
    _X /= _max
    n, d = _X.shape
    rng = np.random.default_rng(random_state)

    centers, _ = kmeans_plusplus(_X, n_clusters=k, random_state=random_state)
    centers = centers.astype(np.float32)

    targets = _targets(n, k, sigma=chain_len_sigma, rng=rng)
    # Build the transport model once
    model, Xvar = _build_transport_lp(n, k, targets, method=1)

    last_shift = np.inf
    start = None  # warm-start matrix for Gurobi

    for it in range(max_iter):
        # 1) Compute cost matrix (squared distances), as float64 for the solver
        cost = _squared_distances(_X, centers).astype(np.float64)

        # 2) Solve transportation LP
        chains, obj, Xmat = _solve_transport_once(model, Xvar, cost, start_X=start)

        # 3) Update centers
        new_centers = np.zeros_like(centers)
        for j in range(k):
            mask = (chains == j)
            # mask must have exactly targets[j] points due to constraints
            new_centers[j] = _X[mask].mean(axis=0) if mask.any() else centers[j]

        shift = float(np.linalg.norm(new_centers - centers))
        centers = new_centers

        # Prepare warm start for next iteration
        start = np.zeros((n, k), dtype=np.float64)
        start[np.arange(n), chains] = 1.0

        if abs(last_shift - shift) < tol:
            break
        last_shift = shift

    # adjust centers
    centers *= _max
    centers += _min
    return chains, centers.astype(np.float32)



def _targets(n, k, sigma=0.0, rng=None):
    """
    Return target sizes per cluster that sum to n.
    - If sigma=0: split as evenly as possible (diff ≤ 1).
    - If sigma>0: allow diversity around the balanced split.

    Parameters
    ----------
    n : int
        Total number of items.
    k : int
        Number of clusters.
    sigma : float
        Standard deviation of Gaussian noise controlling diversity.
    rng : np.random.Generator or None
        Random generator (for reproducibility).

    Returns
    -------
    targets : np.ndarray shape (k,)
        Cluster sizes, sum to n.
    """
    if rng is None:
        rng = np.random.default_rng()

    q, r = divmod(n, k)
    base = np.array([q + 1] * r + [q] * (k - r), dtype=float)

    if sigma > 0:
        noise = rng.normal(0, sigma, size=k)
        noisy = base + noise
        noisy = np.clip(noisy, 1, None)  # no cluster smaller than 1
        # rescale to sum n
        noisy *= n / noisy.sum()
        # round to ints while preserving sum
        targets = np.floor(noisy).astype(int)
        diff = n - targets.sum()
        # fix rounding gap by distributing +1 to some entries
        for i in rng.choice(k, size=diff, replace=False):
            targets[i] += 1
    else:
        targets = base.astype(int)

    return targets


def _build_transport_lp(n: int, k: int, targets: np.ndarray, method=1):
    """
    Build a dense transportation LP:
      min sum_{i,j} c_ij * X_ij
      s.t.  sum_j X_ij = 1  for all i
            sum_i X_ij = targets[j] for all j
            0 <= X_ij <= 1
    Returns (model, X_var) where X_var is an (n,k) MVar.
    """
    m = gp.Model()
    m.Params.OutputFlag = 0
    if method is not None:
        m.Params.Method = method   # 1 = dual simplex (great for repeated solves)

    X = m.addMVar((n, k), lb=0.0, ub=1.0, name="X")

    # Row sums: each point assigned exactly once
    m.addConstrs((X[i, :].sum() == 1 for i in range(n)), name="row")
    # Column sums: exact cluster sizes (integers summing to n)
    m.addConstrs((X[:, j].sum() == int(targets[j]) for j in range(k)), name="col")

    # Dummy objective for now; we set coefficients later each iteration
    m.setObjective(0.0, GRB.MINIMIZE)
    m.update()
    return m, X

def _solve_transport_once(model, X, cost, start_X=None):
    """
    Update objective with new costs and solve.
    cost: (n,k) ndarray of squared distances (float)
    start_X: optional (n,k) ndarray 0/1 warm-start
    Returns labels (n,), objective value, X_matrix (n,k).
    """
    n, k = cost.shape
    # Set objective coefficients without rebuilding the model
    # Easiest: assign Obj on the MVar view (works in recent gurobipy)
    X.Obj = cost
    # Optional warm start
    if start_X is not None:
        X.Start = start_X

    model.optimize()

    # Extract solution
    Xv = X.X  # (n,k)
    labels = Xv.argmax(axis=1).astype(int)
    return labels, model.ObjVal, Xv

def _squared_distances(X, C):
    # GEMM trick: ||x||^2 + ||c||^2 - 2 x·c
    X2 = (X * X).sum(axis=1, keepdims=True)         # (n,1)
    C2 = (C * C).sum(axis=1, keepdims=True).T       # (1,k)
    XC = X @ C.T                                    # (n,k)
    D = X2 + C2 - 2.0 * XC
    # Numerical hygiene: avoid tiny negatives from roundoff
    np.maximum(D, 0, out=D)
    return D

