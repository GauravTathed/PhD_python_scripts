import numpy as np
from numpy.linalg import norm
from collections import deque
import math
import itertools

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

def compute_phase_gate_from_V(V):
    diag = np.diag(V)
    P = np.diag(np.exp(-1j * np.angle(diag)))
    return P

def inverse_single_pulse(U, tol=1e-8):
    U = np.asarray(U, dtype=complex)
    dim = U.shape[0]
    cand = []
    for p in range(dim):
        for q in range(p + 1, dim):
            if abs(U[p, q]) > tol or abs(U[q, p]) > tol:
                cand.append((p, q))
    if len(cand) != 1:
        return (0, 0), 0.0, 0.0
    i, j = cand[0]
    mask = np.ones_like(U, dtype=bool)
    mask[[i, j], :] = False
    mask[:, [i, j]] = False
    if norm(U[mask] - np.eye(dim, dtype=complex)[mask]) > tol:
        raise ValueError("Extra couplings detected")
    a = U[i, i]
    b = U[i, j]
    if not (np.allclose(U[j, j], np.conjugate(a), atol=tol) and np.allclose(U[j, i], -np.conjugate(b), atol=tol)):
        raise ValueError("Not a single two-level unitary of the expected form")
    c = a
    s = b
    c_abs = np.clip(np.abs(c), 0.0, 1.0)
    s_abs = np.clip(np.abs(s), 0.0, 1.0)
    if c_abs == 0.0 and s_abs == 0.0:
        return (i, j), 0.0, 0.0
    theta = 2.0 * np.arctan2(s_abs, c_abs)
    phi = (np.angle(s) - np.angle(c)) % (2 * np.pi)
    return (i, j), float(theta / np.pi), float(phi)

def haar_unitary(n, rng=None):
    rng = np.random.default_rng(rng)
    X = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(X)
    d = np.diag(R)
    D = d / np.abs(d)
    return Q * D

def haar_su(n, rng=None):
    U = haar_unitary(n, rng)
    phi = np.angle(np.linalg.det(U)) / n
    return U * np.exp(-1j * phi)

def scheme_permutations(scheme):
    pair_lists = [pairs for (_, pairs) in scheme]
    perms_per_col = [list(itertools.permutations(pl)) if len(pl) > 1 else [pl] for pl in pair_lists]
    idx_lists = [list(range(len(pl))) for pl in perms_per_col]
    for combo in itertools.product(*idx_lists):
        new_scheme = []
        for (col, _), idx, perms in zip(scheme, combo, perms_per_col):
            new_scheme.append((col, list(perms[idx])))
        yield new_scheme

def decompose_and_reconstruct(U, scheme, tol=1e-10):
    rotation_mats, V = taqr_eliminate(U, scheme, tol)
    P = compute_phase_gate_from_V(V)
    U_rec = P.conj().T
    for G in reversed(rotation_mats):
        U_rec = G.conj().T @ U_rec
    diff = norm(U - U_rec)
    return rotation_mats, P, U_rec, diff

def print_sequence(rotation_mats, P):
    coupling = []
    theta = []
    phases = []
    for k, G in enumerate(rotation_mats):
        (i, j), frac, phi = inverse_single_pulse(G)
        coupling.append((i, j))
        theta.append(frac)
        phases.append(phi)
        print(f"{k:2d}: coupling ({i},{j}), theta={frac:.6f}, phi={phi/np.pi:.6f}")
    phases = np.angle(np.diag(P))
    print("Diagonal phases (per level, /pi):", phases / np.pi)
    return coupling, theta, phases

def demo_perm(dim = 4, rng = 2):
    adj = build_star_adj(dim, center=0)
    scheme = build_taqr_scheme(adj)
    # print("Base TAQR scheme:")
    # print(scheme)
    U = haar_su(dim, rng)
    # QFT = np.array([[1,1,1,1],[1,1j,-1,-1j],[1,-1,1,-1],[1,-1j,-1,1j]], dtype=complex)/2
    # Hadamard = (1/np.sqrt(2))*np.array([[1,1],[1,-1]], dtype=complex)
    # U = np.kron(Hadamard, Hadamard)
    # U = np.kron(U, Hadamard) 
    # U = QFT
    # print("Original U:")
    # print(np.round(U, 3))
    base_rots, base_P, base_Urec, base_diff = decompose_and_reconstruct(U, scheme)
    print("\nBase scheme reconstruction:")
    print("Reconstructed U:")
    print(np.round(base_Urec, 3))
    print("||U - U_rec|| =", base_diff)
    print("Number of rotations:", len(base_rots))
    print("Base sequence:")
    print_sequence(base_rots, base_P)

    print("\nSearching permutations for best scheme...")
    best_scheme = None
    best_rots = None
    best_P = None
    best_diff = None
    best_num = None
    best_theta_sum = None
    # for sch in scheme_permutations(scheme):
    for sch in [scheme]:
        print(sch)
        rots, P, Urec, diff = decompose_and_reconstruct(U, sch)
        if diff >= 1e-8:
            continue
        num = len(rots)
        theta_sum = 0.0
        for G in rots:
            # print(G)
            (_, _), frac, _ = inverse_single_pulse(G)
            theta_sum += frac
        print(theta_sum)
        if best_scheme is None:
            best_scheme = sch
            best_rots = rots
            best_P = P
            best_diff = diff
            best_num = num
            best_theta_sum = theta_sum
        else:
            if num < best_num or (num == best_num and theta_sum < best_theta_sum):
                best_scheme = sch
                best_rots = rots
                best_P = P
                best_diff = diff
                best_num = num
                best_theta_sum = theta_sum

    if best_scheme is not None:
        print("\nBest scheme:")
        print("Scheme:", best_scheme)
        print("||U - U_rec|| =", best_diff)
        print("Number of rotations:", best_num)
        print("Sum theta/pi:", best_theta_sum)
        print("Best sequence:")
        coupling, theta, phase = print_sequence(best_rots, best_P)
        print("Couplings:", coupling)
        print("Thetas (in units of pi):", theta)
        print("Phases (radians):", phase)
    if best_scheme is None:
        print("No valid permutation found.")
    return coupling, theta, phase

print(demo_perm(4,9))