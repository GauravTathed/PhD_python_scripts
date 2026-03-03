import numpy as np
from numpy.linalg import norm
from collections import deque
import math
from scipy.linalg import expm


def coupling_operator_with_phase(i, j, dim, phi):
    op = np.zeros((dim, dim), dtype=complex)
    op[i, j] = np.exp(+1j * phi)
    op[j, i] = np.exp(-1j * phi)
    return op


def evolver_single_pulse(levels, frac, phi, dim, Omega=1.0):
    i, j = levels
    H = 0.5 * Omega * coupling_operator_with_phase(i, j, dim, phi)
    t = (np.pi * frac) / Omega
    return expm(-1j * H * t)


def two_level_left(m, p, q, c, s):
    U = np.eye(m, dtype=complex)
    if p == q:
        return U
    U[p, p] = c
    U[p, q] = s
    U[q, p] = -np.conjugate(s)
    U[q, q] = np.conjugate(c)
    return U


def zero_pair(m, p, q, a, b, tol=1e-12):
    r = np.sqrt(np.abs(a) ** 2 + np.abs(b) ** 2)
    if r < tol:
        return np.eye(m, dtype=complex)
    c = np.conjugate(a) / r
    s = np.conjugate(b) / r
    return two_level_left(m, p, q, c, s)


def bfs(adj, start, active):
    dist = {v: math.inf for v in active}
    parent = {v: None for v in active}
    dist[start] = 0
    dq = deque([start])
    while dq:
        v = dq.popleft()
        for w in adj[v]:
            if w not in active:
                continue
            if dist[w] is math.inf:
                dist[w] = dist[v] + 1
                parent[w] = v
                dq.append(w)
    return dist, parent


def is_connected(adj, active):
    active = set(active)
    if not active:
        return True
    start = next(iter(active))
    dist, _ = bfs(adj, start, active)
    return all(dist[v] is not math.inf for v in active)


def is_non_articulation(adj, active, v):
    active = set(active)
    if v not in active:
        return False
    remaining = active - {v}
    return is_connected(adj, remaining)


def build_taqr_scheme(adj):
    n = len(adj)
    active = set(range(n))
    scheme = []
    while len(active) > 1:
        candidates = [v for v in active if is_non_articulation(adj, active, v)]
        if not candidates:
            raise RuntimeError("No non-articulation vertex found")
        r = max(candidates)
        dist, parent = bfs(adj, r, active)
        L = max(dist[v] for v in active if dist[v] is not math.inf)
        pairs = []
        for level in range(L, 0, -1):
            layer = [v for v in active if dist[v] == level]
            for z in sorted(layer):
                p = parent[z]
                if p is None:
                    continue
                pairs.append((z, p))
        scheme.append((r, pairs))
        active.remove(r)
    remaining = next(iter(active))
    scheme.append((remaining, []))
    return scheme


def build_star_adj(dim, center=0):
    adj = {i: set() for i in range(dim)}
    for j in range(dim):
        if j == center:
            continue
        adj[center].add(j)
        adj[j].add(center)
    return adj


def taqr_eliminate(U, scheme, tol=1e-10):
    U = np.array(U, dtype=complex)
    d = U.shape[0]
    V = U.copy()
    rotation_mats = []
    for col, pairs in scheme:
        for z, p in pairs:
            b = V[z, col]
            if abs(b) < tol:
                continue
            a = V[p, col]
            G = zero_pair(d, p, z, a, b, tol)
            V = G @ V
            rotation_mats.append(G)
    return rotation_mats, V


def inverse_single_taqr_rotation_to_evolver(G, tol=1e-10):
    G = np.asarray(G, dtype=complex)
    dim = G.shape[0]

    cand = [(p, q) for p in range(dim) for q in range(p + 1, dim)
            if abs(G[p, q]) > tol or abs(G[q, p]) > tol]
    if len(cand) != 1:
        raise ValueError(f"Expected exactly one coupled pair, found {cand}")
    i, j = cand[0]

    mask = np.ones_like(G, dtype=bool)
    mask[[i, j], :] = False
    mask[:, [i, j]] = False
    if norm(G[mask] - np.eye(dim, dtype=complex)[mask]) > 100 * tol:
        raise ValueError("Extra couplings detected in rotation")

    c = G[i, i]
    s = G[i, j]

    if not (np.allclose(G[j, j], np.conjugate(c), atol=100 * tol) and
            np.allclose(G[j, i], -np.conjugate(s), atol=100 * tol)):
        raise ValueError("Not TAQR two_level_left form [[c,s],[-s*,c*]]")

    gamma = np.angle(c)

    c_abs = np.clip(np.abs(c), 0.0, 1.0)
    s_abs = np.clip(np.abs(s), 0.0, 1.0)

    theta = 2.0 * np.arctan2(s_abs, c_abs)         # in [0, pi]
    frac = float(theta / np.pi)

    phi = float((np.angle(s) - gamma + np.pi / 2) % (2 * np.pi))

    return (i, j), frac, phi, float(gamma)


def unitary_to_pulses(U, center=0, tol=1e-10, use_virtual_z_shortening=False):
    """
    Input: unitary U (NxN)
    Output:
      couplings: list[(i,j)]
      fractions: list[theta/pi]
      phases  : list[phi_drive] in radians (feed to your evolver)
      frame   : final per-level virtual-Z frame (radians)
    """
    U = np.asarray(U, dtype=complex)
    if U.ndim != 2 or U.shape[0] != U.shape[1]:
        raise ValueError("U must be square")
    dim = U.shape[0]

    adj = build_star_adj(dim, center=center)
    scheme = build_taqr_scheme(adj)

    rotation_mats, _V = taqr_eliminate(U, scheme, tol=tol)

    couplings = []
    fractions = []
    phases = []
    frame = np.zeros(dim, dtype=float)

    for G in rotation_mats:
        (i, j), frac, phi_phys, gamma = inverse_single_taqr_rotation_to_evolver(G, tol=tol)

        frame[i] = (frame[i] + gamma) % (2 * np.pi)
        frame[j] = (frame[j] - gamma) % (2 * np.pi)

        if use_virtual_z_shortening and frac > 1.0 + 1e-12:
            frac = 2.0 - frac
            phi_phys = (phi_phys - np.pi) % (2 * np.pi)
            frame[i] = (frame[i] + np.pi) % (2 * np.pi)
            frame[j] = (frame[j] + np.pi) % (2 * np.pi)

        phi_drive = (phi_phys - (frame[i] - frame[j])) % (2 * np.pi)

        couplings.append((i, j))
        fractions.append(float(frac))
        phases.append(float(phi_drive))

    return couplings, fractions, phases, frame