
import math
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment, minimize_scalar

# =============================================================================
# 137Ba+ constants and manifold definitions
# =============================================================================

I_NUCLEAR = 1.5
G_L = 1.0
G_S = 2.002_319
# 137Ba nuclear magnetic moment (in nuclear magnetons).  The bare nuclear
# g factor is mu_I/(I*mu_N).  The nuclear contribution enters the atomic
# Zeeman Hamiltonian with the opposite sign to the electronic term when
# written in Bohr-magneton units, so G_I_BOHR below is negative.
NUCLEAR_MAGNETIC_MOMENT_MUN = 0.937_365
NUCLEAR_G_FACTOR = NUCLEAR_MAGNETIC_MOMENT_MUN / I_NUCLEAR

E_CHARGE = 1.602_176_634e-19
M_E = 9.109_383_7139e-31
M_P = 1.672_621_923_69e-27
HBAR_NUMERIC = 1.0
H_NUMERIC = 2.0 * np.pi
MU_B_OVER_H = E_CHARGE / (2.0 * M_E) / H_NUMERIC   # Hz/T in the notebook's angular-momentum convention
GAUSS_TO_TESLA = 1e-4
G_I_BOHR = -NUCLEAR_G_FACTOR * (M_E / M_P)

# Hyperfine constants are in Hz.
#
# D3/2 and D5/2 use the measured "uncorrected" constants from Lewty et al.
# These are intentionally used as single-manifold effective constants because
# second-order D3/2-D5/2 hyperfine mixing is not included explicitly here.
# P-state A/B values use the Villemoes et al. experimental values.
#
# Level centers (cm^-1) are NIST Ba II values.
MANIFOLDS = {
    "S12": {
        "label": r"6S1/2",
        "I": 1.5, "J": 0.5, "L": 0, "S": 0.5,
        "A": 4_018.870_833_85e6, "B": 0.0, "C": 0.0,
        "gJ": 2.002_491_92,
        "center_cm": 0.0,
    },
    "P12": {
        "label": r"6P1/2",
        "I": 1.5, "J": 0.5, "L": 1, "S": 0.5,
        "A": 743.7e6, "B": 0.0, "C": 0.0,
        "gJ": None,
        "center_cm": 20_261.561,
    },
    "P32": {
        "label": r"6P3/2",
        "I": 1.5, "J": 1.5, "L": 1, "S": 0.5,
        "A": 127.2e6, "B": 92.5e6, "C": 0.0,
        "gJ": None,
        "center_cm": 21_952.404,
    },
    "D32": {
        "label": r"5D3/2",
        "I": 1.5, "J": 1.5, "L": 2, "S": 0.5,
        "A": 189_730_524.90, "B": 44_538_793.6, "C": 32.465,
        "gJ": 0.799_327_8,
        "center_cm": 4_873.852,
    },
    "D52": {
        "label": r"5D5/2",
        "I": 1.5, "J": 2.5, "L": 2, "S": 0.5,
        "A": -12_029_724.1, "B": 59_519_566.2, "C": -41.73,
        "gJ": 1.200_367_39,
        "center_cm": 5_674.807,
    },
}

_ALIASES = {
    "S1/2": "S12", "6S1/2": "S12", "S12": "S12",
    "P1/2": "P12", "6P1/2": "P12", "P12": "P12",
    "P3/2": "P32", "6P3/2": "P32", "P32": "P32",
    "D3/2": "D32", "5D3/2": "D32", "D32": "D32",
    "D5/2": "D52", "5D5/2": "D52", "D52": "D52",
}

def _name(manifold):
    key = str(manifold).replace("_", "").replace(" ", "")
    for alias, canonical in _ALIASES.items():
        if key.lower() == alias.replace("_", "").replace(" ", "").lower():
            return canonical
    raise ValueError(f"Unknown manifold {manifold!r}. Use one of {list(MANIFOLDS)}.")

# =============================================================================
# Angular momentum and basis helpers
# =============================================================================

def spin_matrices(S):
    """Return Sx, Sy, Sz, S+, S- in the basis m=-S,...,+S."""
    m_vals = np.arange(-S, S + 1, 1.0)
    dim = len(m_vals)
    Sz = np.diag(m_vals.astype(float))
    Sp = np.zeros((dim, dim), dtype=complex)
    Sm = np.zeros((dim, dim), dtype=complex)
    for k, m in enumerate(m_vals[:-1]):
        coef = np.sqrt(S * (S + 1) - m * (m + 1))
        Sp[k + 1, k] = coef
        Sm[k, k + 1] = coef
    Sx = 0.5 * (Sp + Sm)
    Sy = -0.5j * (Sp - Sm)
    return Sx, Sy, Sz, Sp, Sm

def clebschgordan1(j1, m1, j2, m2, j, m):
    """Clebsch-Gordan coefficient <j1 m1; j2 m2 | j m>."""
    vals = (j1, m1, j2, m2, j, m)
    if any(abs(2*x - round(2*x)) > 1e-10 for x in vals):
        raise ValueError("j and m values must be integer or half-integer.")
    if j1 < 0 or j2 < 0 or j < 0 or abs(m1) > j1 or abs(m2) > j2 or abs(m) > j:
        return 0.0
    if abs((m1 + m2) - m) > 1e-10 or j < abs(j1-j2) or j > j1+j2:
        return 0.0

    k_min = max(0, int(round(j2-j-m1)), int(round(j1-j+m2)))
    k_max = min(int(round(j1+j2-j)), int(round(j1-m1)), int(round(j2+m2)))
    if k_min > k_max:
        return 0.0

    fact = lambda x: math.factorial(int(round(x)))
    pref = math.sqrt(
        (2*j+1)
        * fact(j+j1-j2) * fact(j+j2-j1) * fact(j1+j2-j)
        / fact(j1+j2+j+1)
    )
    pref *= math.sqrt(
        fact(j+m) * fact(j-m)
        * fact(j1+m1) * fact(j1-m1)
        * fact(j2+m2) * fact(j2-m2)
    )
    s = 0.0
    for k in range(k_min, k_max + 1):
        denom = (
            fact(k) * fact(j1+j2-j-k) * fact(j1-m1-k)
            * fact(j2+m2-k) * fact(j-j2+m1+k)
            * fact(j-j1-m2+k)
        )
        s += ((-1)**k) / denom
    return pref * s

def uncoupled_basis(I, J):
    """Canonical Hamiltonian basis: |mI,mJ>, mI outer loop, mJ inner loop."""
    return [
        (float(mI), float(mJ))
        for mI in np.arange(-I, I + 1, 1.0)
        for mJ in np.arange(-J, J + 1, 1.0)
    ]

def fmf_labels(I, J):
    """Canonical adiabatic label order: F ascending, then mF ascending."""
    return [
        (int(F), int(mF))
        for F in range(int(abs(I-J)), int(I+J)+1)
        for mF in range(-F, F+1)
    ]

@lru_cache(maxsize=None)
def _fmf_basis_cached(I, J):
    basis = uncoupled_basis(I, J)
    labels = fmf_labels(I, J)
    T = np.zeros((len(basis), len(labels)), dtype=complex)
    for row, (mI, mJ) in enumerate(basis):
        for col, (F, mF) in enumerate(labels):
            T[row, col] = clebschgordan1(I, mI, J, mJ, F, mF)
    return T, tuple(labels)

def make_FmF_basis(I, J):
    """Return T and labels, where T[:,k] is the exact |F,mF> state."""
    T, labels = _fmf_basis_cached(float(I), float(J))
    return T.copy(), list(labels)

# Legacy name, now with unambiguous physical (I,J) argument order.
def transformToFMFBasis(I, J):
    T, labels = make_FmF_basis(I, J)
    rows = [f"mI={mI:+.1f},mJ={mJ:+.1f}" for mI, mJ in uncoupled_basis(I, J)]
    cols = [f"F={F},mF={mF:+d}" for F, mF in labels]
    return T, pd.DataFrame(T, index=rows, columns=cols)

# =============================================================================
# Hamiltonians
# =============================================================================

def lande_g(J, L, S):
    return (
        G_L * (J*(J+1) - S*(S+1) + L*(L+1))
        + G_S * (J*(J+1) + S*(S+1) - L*(L+1))
    ) / (2 * J * (J+1))

def hamiltonian(manifold, B_gauss=0.0):
    """Hyperfine + Zeeman Hamiltonian in Hz in the canonical |mI,mJ> basis."""
    m = MANIFOLDS[_name(manifold)]
    I, J, L, S = m["I"], m["J"], m["L"], m["S"]
    A, Bq, C = m["A"], m["B"], m["C"]

    Ii = np.eye(int(2*I + 1))
    Ij = np.eye(int(2*J + 1))
    Ix, Iy, Iz, _, _ = spin_matrices(I)
    Jx, Jy, Jz, _, _ = spin_matrices(J)

    IJ = np.kron(Ix, Jx) + np.kron(Iy, Jy) + np.kron(Iz, Jz)
    ident = np.eye(IJ.shape[0], dtype=complex)
    IJ2 = IJ @ IJ

    H = A * IJ

    if abs(Bq) > 0 and I >= 1 and J >= 1:
        H = H + Bq * (
            3*IJ2 + 1.5*IJ - I*J*(I+1)*(J+1)*ident
        ) / (2*I*J*(2*I-1)*(2*J-1))

    if abs(C) > 0 and I >= 1.5 and J >= 1.5:
        IJ3 = IJ2 @ IJ
        num = (
            10*IJ3 + 20*IJ2
            + 2*(I*(I+1) + J*(J+1) + 3 - 3*I*(I+1)*J*(J+1))*IJ
            - 5*I*(I+1)*J*(J+1)*ident
        )
        den = I*(I-1)*(2*I-1) * J*(J-1)*(2*J-1)
        H = H + C * num / den

    gJ = m["gJ"] if m.get("gJ") is not None else lande_g(J, L, S)
    B_tesla = float(B_gauss) * GAUSS_TO_TESLA
    H = H + MU_B_OVER_H * B_tesla * (
        gJ * np.kron(Ii, Jz) + G_I_BOHR * np.kron(Iz, Ij)
    )

    return np.real_if_close(H)

def Hamiltonian_S12(B_gauss=0.0): return hamiltonian("S12", B_gauss)
def Hamiltonian_P12(B_gauss=0.0): return hamiltonian("P12", B_gauss)
def Hamiltonian_P32(B_gauss=0.0): return hamiltonian("P32", B_gauss)
def Hamiltonian_D32(B_gauss=0.0): return hamiltonian("D32", B_gauss)
def Hamiltonian_D52(B_gauss=0.0): return hamiltonian("D52", B_gauss)

def solveAndSort(H):
    """Compatibility helper: eigenvectors plus diagonal energy matrix."""
    E, V = np.linalg.eigh(H)
    return V, np.diag(E)

# =============================================================================
# Adiabatic state tracking and labels
# =============================================================================

_EIGEN_CACHE = {}

def clear_eigensystem_cache():
    _EIGEN_CACHE.clear()

def _phase_align(V_old, V_new):
    diag = np.diag(V_old.conj().T @ V_new)
    phases = np.ones_like(diag, dtype=complex)
    nz = np.abs(diag) > 0
    phases[nz] = np.exp(-1j * np.angle(diag[nz]))
    return V_new * phases[np.newaxis, :]

def eigensystem(manifold, B_gauss=0.0, max_step_gauss=0.01, min_overlap=0.97):
    """
    Return a field-dressed eigensystem in fixed adiabatic B->0 |F,mF> order.

    D5/2 labels are *not* reassigned from instantaneous dominant F character.
    They remain attached to the continuous state originating at the zero-field
    |F,mF> state, which is essential because F=3 and F=4 mix strongly at low B.
    """
    name = _name(manifold)
    B_target = float(B_gauss)
    key = (name, round(B_target, 12), float(max_step_gauss), float(min_overlap))
    if key in _EIGEN_CACHE:
        r = _EIGEN_CACHE[key]
        return {
            k: (v.copy() if isinstance(v, np.ndarray) else v.copy() if isinstance(v, pd.DataFrame) else list(v) if isinstance(v, list) else v)
            for k, v in r.items()
        }

    p = MANIFOLDS[name]
    I, J = p["I"], p["J"]
    T, labels = make_FmF_basis(I, J)

    if B_target == 0.0:
        V = T.copy()
    else:
        direction = 1.0 if B_target > 0 else -1.0
        current_B = 0.0
        V = T.copy()

        # Begin very close to zero, then adapt the step. Small initial step is
        # especially important for the near-degenerate D5/2 F=3/F=4 manifold.
        step = min(abs(B_target), min(1e-3, max_step_gauss))
        step = max(step, min(abs(B_target), 1e-8))

        while abs(B_target - current_B) > 1e-14:
            remaining = abs(B_target - current_B)
            trial_step = min(step, remaining)
            next_B = current_B + direction * trial_step

            E_raw, V_raw = np.linalg.eigh(hamiltonian(name, next_B))
            overlap = np.abs(V.conj().T @ V_raw)**2
            old_rows, new_cols = linear_sum_assignment(-overlap)

            col_for_old = np.empty(V.shape[1], dtype=int)
            col_for_old[old_rows] = new_cols
            V_trial = V_raw[:, col_for_old]
            assigned = overlap[np.arange(V.shape[1]), col_for_old]
            quality = float(np.min(assigned))

            if quality < min_overlap and trial_step > 1e-8:
                step = trial_step / 2.0
                continue

            V_trial = _phase_align(V, V_trial)
            V = V_trial
            current_B = next_B

            # Grow where the branches are smooth, shrink automatically where
            # mixing is rapid.
            if quality > 0.9995:
                step = min(max_step_gauss, trial_step * 1.5)
            elif quality < 0.99:
                step = max(1e-8, trial_step * 0.7)
            else:
                step = min(max_step_gauss, trial_step)

    H = hamiltonian(name, B_target)
    energies = np.real(np.diag(V.conj().T @ H @ V))

    # Full instantaneous F,mF decomposition.
    C_fmf = T.conj().T @ V
    P_fmf = np.abs(C_fmf)**2

    # mF is exact for a z-directed static field; use it as an explicit check.
    Ii = np.eye(int(2*I+1))
    Ij = np.eye(int(2*J+1))
    _, _, Iz, _, _ = spin_matrices(I)
    _, _, Jz, _, _ = spin_matrices(J)
    Fz = np.kron(Iz, Ij) + np.kron(Ii, Jz)
    mF_expect = np.real(np.diag(V.conj().T @ Fz @ V))

    rows = []
    for idx, (F0, mF0) in enumerate(labels):
        dom_idx = int(np.argmax(P_fmf[:, idx]))
        Fdom, mFdom = labels[dom_idx]
        rows.append({
            "state_index": idx,
            "label": f"F={F0}, mF={mF0:+d}",
            "F": F0,
            "mF": mF0,
            "energy_Hz": energies[idx],
            "energy_MHz": energies[idx] / 1e6,
            "origin_weight": P_fmf[idx, idx],
            "dominant_F": Fdom,
            "dominant_mF": mFdom,
            "dominant_weight": P_fmf[dom_idx, idx],
            "mF_expectation": mF_expect[idx],
        })
    table = pd.DataFrame(rows)

    result = {
        "manifold": name,
        "B_gauss": B_target,
        "energies_Hz": energies,
        "vectors": V,
        "fmf_coefficients": C_fmf,
        "labels": labels,
        "table": table,
    }
    _EIGEN_CACHE[key] = result
    return {
        k: (v.copy() if isinstance(v, np.ndarray) else v.copy() if isinstance(v, pd.DataFrame) else list(v) if isinstance(v, list) else v)
        for k, v in result.items()
    }

def states(manifold, B_gauss=0.0):
    """Most convenient state-label/energy view."""
    return eigensystem(manifold, B_gauss)["table"]

def state_index(manifold, F, mF):
    labels = fmf_labels(MANIFOLDS[_name(manifold)]["I"], MANIFOLDS[_name(manifold)]["J"])
    try:
        return labels.index((int(F), int(mF)))
    except ValueError as exc:
        raise ValueError(f"{_name(manifold)} has no state F={F}, mF={mF}.") from exc

def state_energy(manifold, F, mF, B_gauss=0.0, unit="MHz"):
    r = eigensystem(manifold, B_gauss)
    E = r["energies_Hz"][state_index(manifold, F, mF)]
    unit_l = unit.lower()
    if unit_l == "hz": return float(E)
    if unit_l == "khz": return float(E/1e3)
    if unit_l == "mhz": return float(E/1e6)
    if unit_l == "ghz": return float(E/1e9)
    raise ValueError("unit must be Hz, kHz, MHz, or GHz.")

# Legacy solver names, but now energies and vectors are guaranteed paired and
# all solvers use the same canonical adiabatic label order.
def _legacy_solver(manifold, B_gauss, diagonal_energy=False):
    r = eigensystem(manifold, B_gauss)
    E = np.diag(r["energies_Hz"]) if diagonal_energy else r["energies_Hz"].copy()
    return E, r["vectors"], r["fmf_coefficients"]

def S12_Energy_at(B_gauss): return _legacy_solver("S12", B_gauss, diagonal_energy=True)
def P12_Energy_at(B_gauss):
    E, V, C = _legacy_solver("P12", B_gauss)
    return E, V, C, np.eye(len(E))
def P32_Energy_at(B_gauss):
    E, V, C = _legacy_solver("P32", B_gauss)
    return E, V, C, np.eye(len(E))
def D32_Energy_at(B_gauss):
    E, V, C = _legacy_solver("D32", B_gauss)
    return E, V, C, np.eye(len(E))
def Energy_at(B_gauss):
    E, V, C = _legacy_solver("D52", B_gauss)
    return E, V, C, np.eye(len(E))
def D52_Energy_at(B_gauss): return Energy_at(B_gauss)

# =============================================================================
# Polarization helpers
# =============================================================================

def elliptical_polarization_from_geometry(
    beam_angle_deg=45.0,
    polarization_angle_deg=0.0,
    ellipticity_angle_deg=0.0,
):
    """
    Electric-field polarization in the physical {sigma-, pi, sigma+} basis.

    Quantization axis: +z.
    Beam: x-z plane, `beam_angle_deg` from +z.
    chi=+45 deg for a beam along +z produces pure sigma+ in this convention.
    chi=-45 deg produces pure sigma-.
    """
    theta = np.deg2rad(beam_angle_deg)
    psi = np.deg2rad(polarization_angle_deg)
    chi = np.deg2rad(ellipticity_angle_deg)

    theta_hat = np.array([np.cos(theta), 0.0, -np.sin(theta)], dtype=complex)
    phi_hat = np.array([0.0, 1.0, 0.0], dtype=complex)
    major = np.cos(psi)*theta_hat + np.sin(psi)*phi_hat
    minor = -np.sin(psi)*theta_hat + np.cos(psi)*phi_hat
    eps = np.cos(chi)*major + 1j*np.sin(chi)*minor
    Ex, Ey, Ez = eps

    # Spherical components of E:
    # E_(+1)=-(Ex+iEy)/sqrt(2) is the physical sigma- field component.
    # E_(-1)= (Ex-iEy)/sqrt(2) is the physical sigma+ field component.
    sigma_minus = -(Ex + 1j*Ey)/np.sqrt(2.0)
    pi = Ez
    sigma_plus = (Ex - 1j*Ey)/np.sqrt(2.0)

    norm = np.sqrt(abs(sigma_minus)**2 + abs(pi)**2 + abs(sigma_plus)**2)
    if norm == 0:
        raise ValueError("Polarization vector has zero norm.")
    return {
        "sigma_minus": sigma_minus/norm,
        "pi": pi/norm,
        "sigma_plus": sigma_plus/norm,
    }

def polarization_components(polarization="pi"):
    """
    Normalize a polarization specification.

    Strings: 'pi', 'sigma+', 'sigma-', 'sigma_plus', 'sigma_minus'.
    A dict may contain sigma_minus, pi, sigma_plus complex amplitudes.
    """
    if isinstance(polarization, str):
        key = polarization.strip().lower().replace(" ", "")
        if key in ("pi", "π"):
            pol = {"sigma_minus": 0j, "pi": 1+0j, "sigma_plus": 0j}
        elif key in ("sigma+", "sigma_plus", "sigmaplus", "σ+"):
            pol = {"sigma_minus": 0j, "pi": 0j, "sigma_plus": 1+0j}
        elif key in ("sigma-", "sigma_minus", "sigmaminus", "σ-"):
            pol = {"sigma_minus": 1+0j, "pi": 0j, "sigma_plus": 0j}
        else:
            raise ValueError("polarization string must be 'pi', 'sigma+', or 'sigma-'.")
    elif isinstance(polarization, dict):
        pol = {
            "sigma_minus": complex(polarization.get("sigma_minus", 0.0)),
            "pi": complex(polarization.get("pi", 0.0)),
            "sigma_plus": complex(polarization.get("sigma_plus", 0.0)),
        }
    else:
        raise TypeError("polarization must be a string or polarization dict.")

    norm = np.sqrt(sum(abs(v)**2 for v in pol.values()))
    if norm == 0:
        raise ValueError("Polarization has zero norm.")
    return {k: v/norm for k, v in pol.items()}

def _pol_for_delta_m(pol):
    """Map transition q=Δm to its physical field component."""
    p = polarization_components(pol)
    # H_int = -d.E and A.B = sum_q (-1)^q A_q B_-q.
    # The global minus in H_int is common to every component; the relative
    # spherical-tensor phase is therefore - for q=+/-1 and + for q=0.
    return {-1: -p["sigma_minus"], 0: p["pi"], +1: -p["sigma_plus"]}

# =============================================================================
# Transition operators and strengths
# =============================================================================

def _e1_uncoupled(initial_manifold, final_manifold, polarization="pi", reduced=1.0):
    ni, nf = _name(initial_manifold), _name(final_manifold)
    pi, pf = MANIFOLDS[ni], MANIFOLDS[nf]
    if pi["I"] != pf["I"]:
        raise ValueError("Initial and final manifolds must have the same nuclear spin.")
    I, Ji, Jf = pi["I"], pi["J"], pf["J"]
    bi, bf = uncoupled_basis(I, Ji), uncoupled_basis(I, Jf)
    eps = _pol_for_delta_m(polarization)
    D = np.zeros((len(bf), len(bi)), dtype=complex)
    for rf, (mIf, mJf) in enumerate(bf):
        for ci, (mIi, mJi) in enumerate(bi):
            if not np.isclose(mIf, mIi):
                continue
            q = int(round(mJf - mJi))
            if q not in (-1, 0, +1):
                continue
            cg = clebschgordan1(Ji, mJi, 1, q, Jf, mJf)
            # <Jf mf|T^1_q|Ji mi> =
            # <Ji mi;1 q|Jf mf> <Jf||T||Ji>/sqrt(2Jf+1)
            D[rf, ci] = reduced * cg / np.sqrt(2*Jf + 1) * eps[q]
    return D

def transition_amplitude(
    initial_manifold,
    final_manifold,
    B_gauss=0.0,
    polarization="pi",
    multipole="E1",
    theta_k=45.0,
    theta_p=58.0,
):
    """
    Complex field-dressed transition-amplitude matrix.
    Rows = final states; columns = initial states.
    """
    ni, nf = _name(initial_manifold), _name(final_manifold)
    Vi = eigensystem(ni, B_gauss)["vectors"]
    Vf = eigensystem(nf, B_gauss)["vectors"]

    kind = multipole.upper()
    if kind == "E1":
        O = _e1_uncoupled(ni, nf, polarization=polarization, reduced=1.0)

    elif kind == "E2":
        # Intended for S1/2 <-> D5/2 relative quadrupole strengths.
        Ji, Jf = MANIFOLDS[ni]["J"], MANIFOLDS[nf]["J"]
        I = MANIFOLDS[ni]["I"]
        bi, bf = uncoupled_basis(I, Ji), uncoupled_basis(I, Jf)
        G = quadrupole_geometric_factors(theta_k, theta_p)
        O = np.zeros((len(bf), len(bi)), dtype=complex)
        for rf, (mIf, mJf) in enumerate(bf):
            for ci, (mIi, mJi) in enumerate(bi):
                if not np.isclose(mIf, mIi):
                    continue
                q = int(round(mJf - mJi))
                if q not in (-2, -1, 0, +1, +2):
                    continue
                cg = clebschgordan1(Ji, mJi, 2, q, Jf, mJf)
                O[rf, ci] = cg * G[q+2]
    else:
        raise ValueError("multipole must be 'E1' or 'E2'.")

    return Vf.conj().T @ O @ Vi

def transition_strength(
    initial_manifold,
    final_manifold,
    B_gauss=0.0,
    polarization="pi",
    multipole="E1",
    theta_k=45.0,
    theta_p=58.0,
    normalize=False,
    as_dataframe=False,
):
    """Physical relative strength |amplitude|^2. Rows final, columns initial."""
    A = transition_amplitude(
        initial_manifold, final_manifold, B_gauss,
        polarization=polarization, multipole=multipole,
        theta_k=theta_k, theta_p=theta_p,
    )
    S = np.abs(A)**2
    if normalize and np.max(S) > 0:
        S = S / np.max(S)
    if as_dataframe:
        li = eigensystem(initial_manifold, B_gauss)["labels"]
        lf = eigensystem(final_manifold, B_gauss)["labels"]
        cols = [f"F={F},mF={m:+d}" for F,m in li]
        rows = [f"F={F},mF={m:+d}" for F,m in lf]
        return pd.DataFrame(S, index=rows, columns=cols)
    return S

def quadrupole_geometric_factors(theta_k=45.0, theta_p=58.0):
    """V1 quadrupole geometric factors, ordered q=-2,-1,0,+1,+2."""
    tk, tp = np.deg2rad(theta_k), np.deg2rad(theta_p)
    return np.array([
        0.25 * abs(np.cos(tp)*np.sin(2*tk) - 2j*np.sin(tp)*np.sin(tk)),
        0.5  * abs(1j*np.sin(tp)*np.cos(tk) + np.cos(tp)*np.cos(2*tk)),
        np.sqrt(6)/4 * abs(np.cos(tp)*np.sin(2*tk)),
        0.5  * abs(1j*np.sin(tp)*np.cos(tk) - np.cos(tp)*np.cos(2*tk)),
        0.25 * abs(np.cos(tp)*np.sin(2*tk) - 2j*np.sin(tp)*np.sin(tk)),
    ], dtype=float)

# Legacy/intuitive wrappers. V2 returns |M|^2 (strength), not |M|.
def TransitionStrength(B_gauss, theta_k=45.0, theta_p=58.0, normalize=False):
    """6S1/2 -> 5D5/2 E2 strength matrix (24 x 8)."""
    return transition_strength(
        "S12", "D52", B_gauss, multipole="E2",
        theta_k=theta_k, theta_p=theta_p, normalize=normalize,
    )

def TransitionStrength_S12_P12(B_gauss, polarization="pi", normalize=False):
    return transition_strength("S12", "P12", B_gauss, polarization, "E1", normalize=normalize)

def TransitionStrength_S12_P32(B_gauss, polarization="pi", normalize=False):
    return transition_strength("S12", "P32", B_gauss, polarization, "E1", normalize=normalize)

def TransitionStrength_D52_P32(B_gauss, polarization="pi", normalize=False):
    return transition_strength("D52", "P32", B_gauss, polarization, "E1", normalize=normalize)

def _magnetic_dipole_uncoupled(manifold, polarization="pi"):
    """M1 operator gJ*J + gI*I, projected onto physical sigma/pi components."""
    name = _name(manifold)
    p = MANIFOLDS[name]
    I, J = p["I"], p["J"]
    gJ = p["gJ"] if p.get("gJ") is not None else lande_g(J, p["L"], p["S"])
    Ix, Iy, Iz, _, _ = spin_matrices(I)
    Jx, Jy, Jz, _, _ = spin_matrices(J)
    Ii, Ij = np.eye(int(2*I+1)), np.eye(int(2*J+1))
    Mx = gJ*np.kron(Ii, Jx) + G_I_BOHR*np.kron(Ix, Ij)
    My = gJ*np.kron(Ii, Jy) + G_I_BOHR*np.kron(Iy, Ij)
    Mz = gJ*np.kron(Ii, Jz) + G_I_BOHR*np.kron(Iz, Ij)

    p = polarization_components(polarization)
    # Interaction B dot mu = sum_q (-1)^q B_q mu_-q.
    # Equivalent Cartesian reconstruction from the physical components:
    # physical sigma+ = B_-1, sigma- = B_+1.
    Bp = p["sigma_minus"]  # spherical +1
    B0 = p["pi"]
    Bm = p["sigma_plus"]   # spherical -1
    Bx = (Bm - Bp)/np.sqrt(2.0)
    By = 1j*(Bm + Bp)/np.sqrt(2.0)
    Bz = B0
    return Bx*Mx + By*My + Bz*Mz

def TransitionStrength_D52_D52(B_gauss, polarization="pi", normalize=False):
    """RF magnetic-dipole strengths within 5D5/2; rows final, columns initial."""
    V = eigensystem("D52", B_gauss)["vectors"]
    O = _magnetic_dipole_uncoupled("D52", polarization)
    A = V.conj().T @ O @ V
    S = np.abs(A)**2
    if normalize and np.max(S)>0:
        S /= np.max(S)
    return S

# =============================================================================
# S1/2 <-> D5/2 optical transition frequencies (legacy table made label-safe)
# =============================================================================

D52_FREQ_ROW_LABELS = (
    [(1,m) for m in range(1,-2,-1)]
    + [(2,m) for m in range(2,-3,-1)]
    + [(3,m) for m in range(3,-4,-1)]
    + [(4,m) for m in range(4,-5,-1)]
)
S12_FREQ_COL_LABELS = [(2,m) for m in range(-2,3)]

def Calculate_raw_freqs(B_gauss):
    d = eigensystem("D52", B_gauss)
    s = eigensystem("S12", B_gauss)
    out = np.empty((24,5), dtype=float)
    for r, (Fd, md) in enumerate(D52_FREQ_ROW_LABELS):
        Ed = d["energies_Hz"][state_index("D52", Fd, md)]
        for c, (Fs, ms) in enumerate(S12_FREQ_COL_LABELS):
            Es = s["energies_Hz"][state_index("S12", Fs, ms)]
            out[r,c] = abs(Ed-Es)/1e6
    return out

def generate_frequencies_at(B_gauss, f0=546.1206708):
    """
    Legacy calibrated 1762-nm relative-frequency table in MHz.
    Reference is D5/2(F=2,mF=0) <-> S1/2(F=2,mF=0), set to f0.
    Forbidden |ΔmF|>2 entries are NaN.
    """
    raw = Calculate_raw_freqs(B_gauss)
    ref_row = D52_FREQ_ROW_LABELS.index((2,0))
    ref_col = S12_FREQ_COL_LABELS.index((2,0))
    out = raw - (raw[ref_row, ref_col] - f0)
    for r, (_, md) in enumerate(D52_FREQ_ROW_LABELS):
        for c, (_, ms) in enumerate(S12_FREQ_COL_LABELS):
            if abs(md-ms) > 2:
                out[r,c] = np.nan
    return out

def fit_B_for_transition(
    f0, f1, row=21, col=2, B_min=0.0, B_max=10.0, tol=1e-9, max_iter=200
):
    def transition_freq(B):
        return generate_frequencies_at(B, f0)[row, col]
    res = minimize_scalar(
        lambda B: (transition_freq(B)-f1)**2,
        bounds=(B_min, B_max), method="bounded",
        options={"xatol": tol, "maxiter": max_iter},
    )
    return float(res.x), float(transition_freq(res.x))

def generate_magnetic_field_sensitivity(B_gauss=4.2095, delta_B=1e-4):
    """d(nu)/dB in MHz/G using a forward finite difference."""
    return (Calculate_raw_freqs(B_gauss + delta_B) - Calculate_raw_freqs(B_gauss)) / delta_B

# =============================================================================
# Raman calculations
# =============================================================================

C_LIGHT = 299_792_458.0

def _center_frequency_hz(lower, upper):
    cm = MANIFOLDS[_name(upper)]["center_cm"] - MANIFOLDS[_name(lower)]["center_cm"]
    return C_LIGHT * 100.0 * cm

def _e1_matrix_eigenbasis(initial, final, B_gauss, polarization, reduced):
    O = _e1_uncoupled(initial, final, polarization, reduced)
    Vi = eigensystem(initial, B_gauss)["vectors"]
    Vf = eigensystem(final, B_gauss)["vectors"]
    return Vf.conj().T @ O @ Vi

def RamanTransitionStrength_S12(
    B_gauss,
    pol1="pi",
    pol2="pi",
    wavelength_nm=532.0,
    reduced_dipole_P12_au=3.3227,
    reduced_dipole_P32_au=4.7017,
    zero_diagonal=False,
):
    """
    Far-detuned Raman coupling inside 6S1/2 through both 6P1/2 and 6P3/2.
    Rows final S states, columns initial S states.
    """
    pol1, pol2 = polarization_components(pol1), polarization_components(pol2)
    nu_laser = C_LIGHT/(wavelength_nm*1e-9)
    rs = eigensystem("S12", B_gauss)
    E_S = rs["energies_Hz"]

    amplitudes = {}
    for P, red in [("P12", reduced_dipole_P12_au), ("P32", reduced_dipole_P32_au)]:
        rp = eigensystem(P, B_gauss)
        E_P = rp["energies_Hz"]
        nu0 = _center_frequency_hz("S12", P)
        D1 = _e1_matrix_eigenbasis("S12", P, B_gauss, pol1, red)
        D2 = _e1_matrix_eigenbasis("S12", P, B_gauss, pol2, red)
        A = np.zeros((len(E_S), len(E_S)), dtype=complex)
        for i in range(len(E_S)):
            det = nu_laser - (nu0 + E_P - E_S[i])
            if np.any(np.abs(det) < 1.0):
                raise ValueError("Laser too close to an intermediate resonance.")
            A[:,i] = D2.conj().T @ (D1[:,i]/det)
        amplitudes[P] = A

    A_total = amplitudes["P12"] + amplitudes["P32"]
    if zero_diagonal:
        np.fill_diagonal(A_total,0)
        for A in amplitudes.values(): np.fill_diagonal(A,0)
    S = abs(A_total)**2
    rel = S/np.max(S) if np.max(S)>0 else S.copy()
    return {
        "amplitude": A_total,
        "strength": S,
        "relative_strength": rel,
        "P12_amplitude": amplitudes["P12"],
        "P32_amplitude": amplitudes["P32"],
        "ground_energies_hz": E_S.copy(),
        "labels": rs["labels"],
    }

def RamanTransitionStrength_D52(
    B_gauss,
    pol1="pi",
    pol2="pi",
    wavelength_nm=532.0,
    reduced_dipole_D52_P32_au=4.1028,
    zero_diagonal=False,
):
    """Far-detuned Raman coupling inside 5D5/2 through 6P3/2."""
    pol1, pol2 = polarization_components(pol1), polarization_components(pol2)
    nu_laser = C_LIGHT/(wavelength_nm*1e-9)
    rd = eigensystem("D52", B_gauss)
    rp = eigensystem("P32", B_gauss)
    E_D, E_P = rd["energies_Hz"], rp["energies_Hz"]
    nu0 = _center_frequency_hz("D52","P32")
    D1 = _e1_matrix_eigenbasis("D52","P32",B_gauss,pol1,reduced_dipole_D52_P32_au)
    D2 = _e1_matrix_eigenbasis("D52","P32",B_gauss,pol2,reduced_dipole_D52_P32_au)

    A = np.zeros((len(E_D),len(E_D)),dtype=complex)
    detmat = np.zeros((len(E_P),len(E_D)))
    for i in range(len(E_D)):
        det = nu_laser - (nu0 + E_P - E_D[i])
        detmat[:,i]=det
        if np.any(np.abs(det)<1.0):
            raise ValueError("Laser too close to an intermediate resonance.")
        A[:,i]=D2.conj().T @ (D1[:,i]/det)
    if zero_diagonal: np.fill_diagonal(A,0)
    S=abs(A)**2
    rel=S/np.max(S) if np.max(S)>0 else S.copy()
    return {
        "amplitude":A, "strength":S, "relative_strength":rel,
        "P32_amplitude":A.copy(), "ground_energies_hz":E_D.copy(),
        "detuning_hz":detmat, "labels":rd["labels"],
        "laser_frequency_hz":nu_laser,
        "D52_to_P32_center_frequency_hz":nu0,
    }

# =============================================================================
# Validation
# =============================================================================

def run_self_tests(verbose=True):
    """Numerical consistency tests for basis, solvers, labels, and polarization."""
    records = []

    def add(test, detail, error, tolerance):
        records.append({
            "test": test, "detail": detail,
            "error": float(error), "tolerance": float(tolerance),
            "passed": bool(error <= tolerance),
        })

    # Hamiltonians and canonical |F,mF> bases.
    for name, p in MANIFOLDS.items():
        H = hamiltonian(name, 4.2079)
        add(name, "Hamiltonian Hermiticity", np.max(np.abs(H-H.conj().T)), 1e-8)

        T, _ = make_FmF_basis(p["I"], p["J"])
        add(name, "FmF basis unitarity",
            np.max(np.abs(T.conj().T@T - np.eye(T.shape[1]))), 1e-12)

        H0f = T.conj().T @ hamiltonian(name,0.0) @ T
        add(name, "H(B=0) diagonal in FmF",
            np.max(np.abs(H0f-np.diag(np.diag(H0f)))), 2e-6)

    # D3/2 zero-field intervals: measured values reproduced by the effective
    # single-manifold constants used above.
    d32_0 = states("D32",0.0).groupby("F")["energy_Hz"].mean()
    measured_d32 = {
        (0,1): 145_193_549.3,
        (1,2): 334_921_347.13,
        (2,3): 613_730_628.08,
    }
    for (Fa,Fb), expected in measured_d32.items():
        got = abs(float(d32_0[Fa]-d32_0[Fb]))
        add("D32", f"zero-field F={Fa}->F={Fb} interval", abs(got-expected), 0.2)

    # D5/2 zero-field intervals: Lewty et al. 2013 measured values.
    d0 = states("D52",0.0).groupby("F")["energy_Hz"].mean()
    measured = {
        (1,2): 71_675_902.4,
        (2,3): 62_872_301.0,
        (3,4): 503_510.5,
    }
    for (Fa,Fb), expected in measured.items():
        got = float(d0[Fa]-d0[Fb])
        add("D52", f"zero-field F={Fa}->F={Fb} interval", abs(got-expected), 0.2)

    # Field-dressed eigensystems.
    for name in MANIFOLDS:
        for B in (0.01, 0.5, 4.2079, 10.0):
            r=eigensystem(name,B)
            V,E=r["vectors"],r["energies_Hz"]
            H=hamiltonian(name,B)
            add(f"{name}@{B}G","orthonormality",
                np.max(np.abs(V.conj().T@V-np.eye(V.shape[1]))),1e-11)
            add(f"{name}@{B}G","eigen residual (Hz)",
                np.max(np.linalg.norm(H@V - V*E[np.newaxis,:],axis=0)),5e-6)
            add(f"{name}@{B}G","mF identity",
                np.max(np.abs(r["table"]["mF_expectation"].to_numpy()
                              -r["table"]["mF"].to_numpy())),1e-10)

    # Polarization convention.
    sp=elliptical_polarization_from_geometry(0,0,+45)
    sm=elliptical_polarization_from_geometry(0,0,-45)
    add("polarization","chi=+45 at k||+z -> sigma+",
        1-abs(sp["sigma_plus"])**2,1e-12)
    add("polarization","chi=-45 at k||+z -> sigma-",
        1-abs(sm["sigma_minus"])**2,1e-12)

    # Exact E1 ΔmF selection rules for all implemented optical dipoles.
    for initial, final in [("S12","P12"),("S12","P32"),("D52","P32")]:
        li=eigensystem(initial,4.2079)["labels"]
        lf=eigensystem(final,4.2079)["labels"]
        for pol,dm in [("sigma+",1),("pi",0),("sigma-",-1)]:
            S=transition_strength(initial,final,4.2079,pol)
            forbidden=0.0
            for f,(_,mf) in enumerate(lf):
                for i,(_,mi) in enumerate(li):
                    if mf-mi != dm:
                        forbidden=max(forbidden,float(S[f,i]))
            add(f"{initial}->{final}",f"{pol} forbidden strength",forbidden,1e-20)

    # M1 ΔmF selection rule.
    labels=eigensystem("D52",4.2079)["labels"]
    for pol,dm in [("sigma+",1),("pi",0),("sigma-",-1)]:
        S=TransitionStrength_D52_D52(4.2079,pol)
        forbidden=0.0
        for f,(_,mf) in enumerate(labels):
            for i,(_,mi) in enumerate(labels):
                if mf-mi != dm:
                    forbidden=max(forbidden,float(S[f,i]))
        add("D52 M1",f"{pol} forbidden strength",forbidden,1e-20)

    df=pd.DataFrame(records)
    passed=bool(df["passed"].all())
    if verbose:
        print("All tests passed." if passed else "Some tests failed.")
        try:
            display(df)
        except NameError:
            print(df.to_string(index=False))
    return passed,df
