import math
from collections import deque

import numpy as np


def unitary_to_G_rotations(U, center=0, tol=1e-10):
    """
    Input: unitary U (NxN)
    Output:
      rotation_mats: list of TAQR elimination matrices
      V: residual matrix after elimination
    """

    def two_level_left(m, p, q, c, s):
        G = np.eye(m, dtype=complex)
        if p == q:
            return G
        G[p, p] = c
        G[p, q] = s
        G[q, p] = -np.conjugate(s)
        G[q, q] = np.conjugate(c)
        return G

    def zero_pair(m, p, q, a, b):
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
        queue = deque([start])
        while queue:
            v = queue.popleft()
            for w in adj[v]:
                if w not in active or dist[w] != math.inf:
                    continue
                dist[w] = dist[v] + 1
                parent[w] = v
                queue.append(w)
        return dist, parent

    def is_connected(adj, active):
        active = set(active)
        if not active:
            return True
        start = next(iter(active))
        dist, _ = bfs(adj, start, active)
        return all(dist[v] != math.inf for v in active)

    def is_non_articulation(adj, active, v):
        active = set(active)
        if v not in active:
            return False
        return is_connected(adj, active - {v})

    def build_taqr_scheme(adj):
        active = set(adj)
        scheme = []
        while len(active) > 1:
            candidates = [v for v in active if is_non_articulation(adj, active, v)]
            if not candidates:
                raise RuntimeError("No non-articulation vertex found")
            root = max(candidates)
            dist, parent = bfs(adj, root, active)
            depth = max(dist[v] for v in active if dist[v] != math.inf)
            pairs = []
            for level in range(depth, 0, -1):
                layer = [v for v in active if dist[v] == level]
                for z in sorted(layer):
                    p = parent[z]
                    if p is not None:
                        pairs.append((z, p))
            scheme.append((root, pairs))
            active.remove(root)
        scheme.append((next(iter(active)), []))
        return scheme

    def build_star_adj(dim):
        adj = {i: set() for i in range(dim)}
        for j in range(dim):
            if j == center:
                continue
            adj[center].add(j)
            adj[j].add(center)
        return adj

    def taqr_eliminate(U_mat, scheme):
        U_mat = np.array(U_mat, dtype=complex)
        dim = U_mat.shape[0]
        V = U_mat.copy()
        rotation_mats = []
        for col, pairs in scheme:
            for z, p in pairs:
                b = V[z, col]
                if abs(b) < tol:
                    continue
                a = V[p, col]
                G = zero_pair(dim, p, z, a, b)
                V = G @ V
                rotation_mats.append(G)
        return rotation_mats, V

    U = np.asarray(U, dtype=complex)
    if U.ndim != 2 or U.shape[0] != U.shape[1]:
        raise ValueError("U must be square")

    adj = build_star_adj(U.shape[0])
    scheme = build_taqr_scheme(adj)
    return taqr_eliminate(U, scheme)