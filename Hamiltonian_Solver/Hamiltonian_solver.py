import numpy as np
import pandas as pd
import math
from scipy.optimize import minimize_scalar

# =============================================================================
# Helper Functions
# =============================================================================

def spin_matrices(S: float):
    m_vals = np.arange(-S, S + 1)          
    dim = len(m_vals)

    Sz = np.diag(m_vals)

    Sp = np.zeros((dim, dim), dtype=complex)
    Sm = np.zeros((dim, dim), dtype=complex)

    for k, m in enumerate(m_vals[:-1]):   
        coef = np.sqrt(S * (S + 1) - m * (m + 1))
        if coef:                        
            Sp[k + 1, k] = coef
            Sm[k, k + 1] = coef         

    Sx = 0.5 * (Sp + Sm)
    Sy = -0.5j * (Sp - Sm)

    return Sx, Sy, Sz, Sp, Sm

def transformToFMFBasis(I, J):
    dim_I = int(2 * I + 1)
    dim_J = int(2 * J + 1)
    dim_F = dim_I * dim_J
    T = np.zeros((dim_F, dim_F))
    col_names = []
    row_names = []
    idx_FmF = 0
    for F in reversed(range(int(abs(J - I)), int(J + I) + 1)):
        for mF in reversed(np.arange(-F, F + 1)):
            idx_FmF += 1
            col_names.append(f'F{F}_mF{mF}')
            idx_IJmImJ = 0
            for mJ in reversed(np.arange(-J, J + 1)):
                for mI in reversed(np.arange(-I, I + 1)):
                    idx_IJmImJ += 1
                    if idx_FmF == 1:
                        row_names.append(f'mI {round(mJ,1)}, mJ {round(mI,1)}')
                    cg_coeff = clebschgordan1(J, mJ, I, mI, F, mF)
                    T[idx_IJmImJ - 1, idx_FmF - 1] = cg_coeff
    T_table = pd.DataFrame(T, index=row_names, columns=col_names)
    return T, T_table

def clebschgordan1(j1,m1,j2,m2,j,m):
    if j1<0 or j2<0 or j<0 or (2*j1)%1 or (2*j2)%1 or (2*j)%1 or (2*m1)%1 or (2*m2)%1 or (2*m)%1 or abs(m1)>j1 or abs(m2)>j2 or abs(m)>j or j1+m1<0 or j2+m2<0 or j+m<0 or j1+j2+j<0 or (j1+m1)%1 or (j2+m2)%1 or (j+m)%1 or (j1+j2+j)%1:
        raise ValueError
    if m1+m2-m or j<abs(j1-j2) or j>j1+j2:
        return 0.0
    k_min=max(0,j2-j-m1,j1-j+m2)
    k_max=min(j1+j2-j,j1-m1,j2+m2)
    pref=math.sqrt((2*j+1)*math.factorial(int(j+j1-j2))*math.factorial(int(j+j2-j1))*math.factorial(int(j1+j2-j))/math.factorial(int(j1+j2+j+1)))
    pref*=math.sqrt(math.factorial(int(j+m))*math.factorial(int(j-m))*math.factorial(int(j1+m1))*math.factorial(int(j1-m1))*math.factorial(int(j2+m2))*math.factorial(int(j2-m2)))
    s=0.0
    for k in range(int(k_min),int(k_max)+1):
        denom=math.factorial(k)*math.factorial(int(j1+j2-j-k))*math.factorial(int(j1-m1-k))*math.factorial(int(j2+m2-k))*math.factorial(int(j-j2+m1+k))*math.factorial(int(j-j1-m2+k))
        s+=((-1)**k)/denom
    return pref*s

def CG_mat_S12(I, J):
    dim_I = int(2 * I + 1)
    dim_J = int(2 * J + 1)
    dim_F = dim_I * dim_J
    T = np.zeros((dim_F, dim_F))
    col_names = []
    row_names = []
    idx_FmF = 0
    F = 1
    for mF in range(-F, F + 1):
        idx_FmF += 1
        col_names.append(f'F{F}_mF{mF}')
        idx_IJmImJ = 0
        for mJ in np.flip(np.arange(-J, J + 1)):
            for mI in np.flip(np.arange(-I, I + 1)):
                idx_IJmImJ += 1
                cg = clebschgordan1(J, mJ, I, mI, F, mF)
                T[idx_IJmImJ - 1, idx_FmF - 1] = cg
    F = 2
    for mF in np.flip(np.arange(-F, F + 1)):
        idx_FmF += 1
        col_names.append(f'F{F}_mF{mF}')
        idx_IJmImJ = 0
        for mJ in np.flip(np.arange(-J, J + 1)):
            for mI in np.flip(np.arange(-I, I + 1)):
                idx_IJmImJ += 1
                if idx_FmF == 6:
                    row_names.append(f'mI {round(mJ, 1)}, mJ {round(mI, 1)}')
                cg = clebschgordan1(J, mJ, I, mI, F, mF)
                T[idx_IJmImJ - 1, idx_FmF - 1] = cg
    T_table = pd.DataFrame(T, index=row_names, columns=col_names)
    return T, T_table

def solveAndSort(Hamiltonian):
    energyLvlUnsorted, vectorsUnsorted = np.linalg.eigh(Hamiltonian)
    ind = np.argsort(energyLvlUnsorted)
    sortedEigenvectors = vectorsUnsorted[:, ind]
    sorted_energy_lvl_D52 = np.diag(energyLvlUnsorted[ind])
    return sortedEigenvectors, sorted_energy_lvl_D52

# =============================================================================
# Hamiltonians
# =============================================================================

def Hamiltonian_D52(B_e_gauss: float) -> np.ndarray:

    I = 3/2
    J = 5/2
    L = 2
    S = 1/2

    g_L = 1.0
    g_S = 2.002_319

    # Hyperfine constants (Hz) from Lewty 2013
    A_D = -12_029_724.1        # magnetic dipole
    B_Q =   59_519_566.2       # electric quadrupole
    C_O =          -41.73      # electric octupole

    e  = 1.602_176_63e-19     
    m_e = 9.109_383_7e-31     
    m_p = 1.672_621_92e-27     
    mu_B = e / (2 * m_e)       
    h = 2 * np.pi              

    g_J = ( g_L * (J*(J+1) - S*(S+1) + L*(L+1))
          + g_S * (J*(J+1) + S*(S+1) - L*(L+1)) ) / (2 * J * (J+1))
    g_I = 0.624_867 * (m_e / m_p)

    gauss_to_tesla = 1e-4
    B_e = B_e_gauss * gauss_to_tesla

    Ii = np.eye(int(2*I + 1))
    Ij = np.eye(int(2*J + 1))

    Ix, Iy, Iz, Ip, Im = spin_matrices(I)
    Jx, Jy, Jz, Jp, Jm = spin_matrices(J)

    IJ = (np.kron(Ix, Jx) + np.kron(Iy, Jy) + np.kron(Iz, Jz))
    ident = np.eye(IJ.shape[0])

    IJ2 = IJ @ IJ
    IJ3 = IJ2 @ IJ

    H_A = A_D * IJ

    H_B = B_Q * ( 3*IJ2 + 1.5*IJ - I*J*(I+1)*(J+1)*ident ) / ( 2*I*J*(2*I-1)*(2*J-1) )

    H_Z = mu_B * B_e * ( g_J*np.kron(Ii, Jz) + g_I*np.kron(Iz, Ij) ) / h

    num = 10*IJ3 + 20*IJ2 + 2*(I*(I+1) + J*(J+1) + 3 - 3*I*(I+1)*J*(J+1)) * IJ \
          - 5*I*(I+1)*J*(J+1) * ident
    den = I*(I-1)*(2*I-1) * J*(J-1)*(2*J-1)
    H_C = C_O * num / den

    return H_A + H_B + H_C + H_Z

def Hamiltonian_S12(Be):
    I = 3 / 2
    J = 1 / 2
    L = 0
    S = 1 / 2
    gL = 1
    gS = 2.002319
    AD = 4018.8708338e6
    ec = 1.60217663e-19
    me = 9.1093837e-31
    mp = 1.67262192e-27
    hbar = 1.054571817e-34
    h = 6.62607015e-34
    mu_B = ec * hbar / (2 * me)
    gJ = (gL * (J * (J + 1) - S * (S + 1) + L * (L + 1)) + gS * (J * (J + 1) + S * (S + 1) - L * (L + 1))) / (2 * J * (J + 1))
    gI = 0.624867 * (me / mp)
    gausstotesla = 1e-4
    Ii = np.eye(int(2 * I + 1))
    Ij = np.eye(int(2 * J + 1))
    Ix, Iy, Iz, Ip, Im = spin_matrices(I)
    Jx, Jy, Jz, Jp, Jm = spin_matrices(J)
    IJ = np.kron(Ix, Jx) + np.kron(Iy, Jy) + np.kron(Iz, Jz)
    H_md = AD * IJ
    H_z = mu_B * Be * gausstotesla * (gJ * np.kron(Ii, Jz) + gI * np.kron(Iz, Ij)) / h
    return H_md + H_z

def Hamiltonian_P32(B_e_gauss):
    I = 3/2
    J = 3/2
    L = 1
    S = 1/2

    g_L = 1.0
    g_S = 2.002_319

    # Hyper-fine constants (Hz) for 137Ba+ 6P3/2: Villemose (1993)
    A_P = 126_900_000.0      
    B_Q = 95_500_000.0       
    C_O = 0.0                

    e  = 1.602_176_63e-19
    m_e = 9.109_383_7e-31
    m_p = 1.672_621_92e-27
    mu_B = e / (2 * m_e)
    h = 2 * np.pi

    g_J = ( g_L * (J*(J+1) - S*(S+1) + L*(L+1))
          + g_S * (J*(J+1) + S*(S+1) - L*(L+1)) ) / (2 * J * (J+1))
    g_I = 0.624_867 * (m_e / m_p)

    gauss_to_tesla = 1e-4
    B_e = B_e_gauss * gauss_to_tesla

    Ii = np.eye(int(2*I + 1))
    Ij = np.eye(int(2*J + 1))

    Ix, Iy, Iz, Ip, Im = spin_matrices(I)
    Jx, Jy, Jz, Jp, Jm = spin_matrices(J)

    IJ = (np.kron(Ix, Jx) + np.kron(Iy, Jy) + np.kron(Iz, Jz))
    ident = np.eye(IJ.shape[0])

    IJ2 = IJ @ IJ
    IJ3 = IJ2 @ IJ

    H_A = A_P * IJ
    H_B = B_Q * ( 3*IJ2 + 1.5*IJ - I*J*(I+1)*(J+1)*ident ) / ( 2*I*J*(2*I-1)*(2*J-1) )
    H_Z = mu_B * B_e * ( g_J*np.kron(Ii, Jz) + g_I*np.kron(Iz, Ij) ) / h
    H_C = C_O * ( 10*IJ3 + 20*IJ2 + 2*(I*(I+1)+J*(J+1)+3-3*I*(I+1)*J*(J+1))*IJ
                 -5*I*(I+1)*J*(J+1)*ident ) / ( I*(I-1)*(2*I-1) * J*(J-1)*(2*J-1) )

    return H_A + H_B + H_C + H_Z

def Hamiltonian_P12(B_e_gauss):
    """
    Renamed from Hamiltonian_P32 in the original notebook.
    Corresponds to P1/2 Hamiltonian based on the notebook header.
    Note: The code uses J=3/2, which might be a copy-paste error in the original notebook
    if it was intended for P1/2.
    """
    I = 3/2
    J = 3/2
    L = 1
    S = 1/2

    g_L = 1.0
    g_S = 2.002_319

    # Hyper-fine constants (Hz) for 137Ba+ 6P3/2: Villemose (1993)
    A_P = 747.7e6     
    B_Q = 0       
    C_O = 0              

    e  = 1.602_176_63e-19
    m_e = 9.109_383_7e-31
    m_p = 1.672_621_92e-27
    mu_B = e / (2 * m_e)
    h = 2 * np.pi

    g_J = ( g_L * (J*(J+1) - S*(S+1) + L*(L+1))
          + g_S * (J*(J+1) + S*(S+1) - L*(L+1)) ) / (2 * J * (J+1))
    g_I = 0.624_867 * (m_e / m_p)

    gauss_to_tesla = 1e-4
    B_e = B_e_gauss * gauss_to_tesla

    Ii = np.eye(int(2*I + 1))
    Ij = np.eye(int(2*J + 1))

    Ix, Iy, Iz, Ip, Im = spin_matrices(I)
    Jx, Jy, Jz, Jp, Jm = spin_matrices(J)

    IJ = (np.kron(Ix, Jx) + np.kron(Iy, Jy) + np.kron(Iz, Jz))
    ident = np.eye(IJ.shape[0])

    IJ2 = IJ @ IJ
    IJ3 = IJ2 @ IJ

    H_A = A_P * IJ
    H_B = B_Q * ( 3*IJ2 + 1.5*IJ - I*J*(I+1)*(J+1)*ident ) / ( 2*I*J*(2*I-1)*(2*J-1) )
    H_Z = mu_B * B_e * ( g_J*np.kron(Ii, Jz) + g_I*np.kron(Iz, Ij) ) / h
    H_C = C_O * ( 10*IJ3 + 20*IJ2 + 2*(I*(I+1)+J*(J+1)+3-3*I*(I+1)*J*(J+1))*IJ
                 -5*I*(I+1)*J*(J+1)*ident ) / ( I*(I-1)*(2*I-1) * J*(J-1)*(2*J-1) )

    return H_A + H_B + H_C + H_Z

def Hamiltonian_D32(B_e_gauss: float) -> np.ndarray:

    I = 3/2
    J = 3/2
    L = 2
    S = 1/2

    g_L = 1.0
    g_S = 2.002_319

    # Hyperfine constants (Hz) from Lewty 2013
    A_D = 189731494        # magnetic dipole
    B_Q = 44537594       # electric quadrupole
    C_O = 29.533      # electric octupole

    e  = 1.602_176_63e-19     
    m_e = 9.109_383_7e-31     
    m_p = 1.672_621_92e-27     
    mu_B = e / (2 * m_e)       
    h = 2 * np.pi              

    g_J = ( g_L * (J*(J+1) - S*(S+1) + L*(L+1))
          + g_S * (J*(J+1) + S*(S+1) - L*(L+1)) ) / (2 * J * (J+1))
    g_I = 0.624_867 * (m_e / m_p)

    gauss_to_tesla = 1e-4
    B_e = B_e_gauss * gauss_to_tesla

    Ii = np.eye(int(2*I + 1))
    Ij = np.eye(int(2*J + 1))

    Ix, Iy, Iz, Ip, Im = spin_matrices(I)
    Jx, Jy, Jz, Jp, Jm = spin_matrices(J)

    IJ = (np.kron(Ix, Jx) + np.kron(Iy, Jy) + np.kron(Iz, Jz))
    ident = np.eye(IJ.shape[0])

    IJ2 = IJ @ IJ
    IJ3 = IJ2 @ IJ

    H_A = A_D * IJ

    H_B = B_Q * ( 3*IJ2 + 1.5*IJ - I*J*(I+1)*(J+1)*ident ) / ( 2*I*J*(2*I-1)*(2*J-1) )

    H_Z = mu_B * B_e * ( g_J*np.kron(Ii, Jz) + g_I*np.kron(Iz, Ij) ) / h

    num = 10*IJ3 + 20*IJ2 + 2*(I*(I+1) + J*(J+1) + 3 - 3*I*(I+1)*J*(J+1)) * IJ \
          - 5*I*(I+1)*J*(J+1) * ident
    den = I*(I-1)*(2*I-1) * J*(J-1)*(2*J-1)
    H_C = C_O * num / den

    return H_A + H_B + H_C + H_Z

# =============================================================================
# Solvers
# =============================================================================

def Energy_at(stop_Be):
    J = 5 / 2
    I = 3 / 2
    prev_Be = 1e-5
    vectors_low, _ = solveAndSort(Hamiltonian_D52(prev_Be))
    T_D52, _ = transformToFMFBasis(J, I)
    iden = np.eye(24)
    vectors_FmF_low = iden.copy()
    position_check_matrix = iden.copy()
    it = 0
    for Be in np.arange(1e-6, stop_Be + 0.0001 * stop_Be + 1e-12, 0.0001 * stop_Be):
        vectors, energy_mat = solveAndSort(Hamiltonian_D52(Be))
        vectors_high = vectors
        pcm_new = np.abs(np.round(vectors_high.T @ vectors_low))
        position_check_matrix = np.abs(np.round(position_check_matrix @ pcm_new))
        vectors = vectors @ position_check_matrix.T
        energy_vec = position_check_matrix @ np.diag(energy_mat)
        if it == 0:
            vectors_D52_FmF = T_D52.T @ vectors
            iden_first = np.diag(np.round(np.sum(vectors_D52_FmF, axis=0)))
            vectors = vectors @ iden_first
        new_iden = np.round(pcm_new)
        iden = iden @ new_iden
        vectors_D52_FmF = T_D52.T @ vectors
        iden_new = np.round(vectors_D52_FmF.T @ vectors_FmF_low)
        vectors = vectors @ iden_new
        vectors_D52_FmF = T_D52.T @ vectors
        vectors[np.abs(vectors) < 1e-6] = 0
        vectors_D52_FmF[np.abs(vectors_D52_FmF) < 1e-6] = 0
        vectors_low = vectors_high
        vectors_FmF_low = vectors_D52_FmF
        it += 1
    return energy_vec, vectors, vectors_D52_FmF, position_check_matrix

def S12_Energy_at(stop_Be):
    J = 0.5
    I = 1.5
    vectors_low, _ = solveAndSort(Hamiltonian_S12(1e-5))
    iden = np.eye(8)
    vectors_FmF_low = iden.copy()
    position_check_matrix = iden.copy()
    it = 0
    T_S12, _ = CG_mat_S12(J, I)
    step = 0.0001 * stop_Be
    for Be in np.arange(1e-7, stop_Be + step + 1e-12, step):
        vectors_S12, energy_lvl_S12 = solveAndSort(Hamiltonian_S12(Be))
        vectors_high = vectors_S12
        pcm_new = np.abs(np.round(vectors_high.T @ vectors_low))
        position_check_matrix = np.abs(np.round(position_check_matrix @ pcm_new))
        vectors_S12 = vectors_S12 @ position_check_matrix.T
        if it == 0:
            vectors_S12_FmF = T_S12.T @ vectors_S12
            iden_first = np.diag(np.round(np.sum(vectors_S12_FmF, axis=0)))
            vectors_S12 = vectors_S12 @ iden_first
        iden = iden @ np.round(pcm_new)
        vectors_S12_FmF = T_S12.T @ vectors_S12
        vectors_S12 = vectors_S12 @ np.round(vectors_S12_FmF.T @ vectors_FmF_low)
        vectors_S12_FmF = T_S12.T @ vectors_S12
        vectors_low = vectors_high
        vectors_FmF_low = vectors_S12_FmF
        it += 1
    return energy_lvl_S12, vectors_S12, vectors_S12_FmF

def P32_Energy_at(stop_Be):
    J = 3 / 2
    I = 3 / 2
    prev_Be = 1e-5
    vectors_low, _ = solveAndSort(Hamiltonian_P32(prev_Be))
    T_P32, _ = transformToFMFBasis(J, I)
    iden = np.eye(16)
    vectors_FmF_low = iden.copy()
    position_check_matrix = iden.copy()
    it = 0
    for Be in np.arange(1e-6, stop_Be + 0.0001 * stop_Be + 1e-12, 0.0001 * stop_Be):
        vectors, energy_mat = solveAndSort(Hamiltonian_P32(Be))
        vectors_high = vectors
        pcm_new = np.abs(np.round(vectors_high.T @ vectors_low))
        position_check_matrix = np.abs(np.round(position_check_matrix @ pcm_new))
        vectors = vectors @ position_check_matrix.T
        energy_vec = position_check_matrix @ np.diag(energy_mat)
        if it == 0:
            vectors_P32_FmF = T_P32.T @ vectors
            iden_first = np.diag(np.round(np.sum(vectors_P32_FmF, axis=0)))
            vectors = vectors @ iden_first
        new_iden = np.round(pcm_new)
        iden = iden @ new_iden
        vectors_P32_FmF = T_P32.T @ vectors
        iden_new = np.round(vectors_P32_FmF.T @ vectors_FmF_low)
        vectors = vectors @ iden_new
        vectors_P32_FmF = T_P32.T @ vectors
        vectors[np.abs(vectors) < 1e-6] = 0
        vectors_P32_FmF[np.abs(vectors_P32_FmF) < 1e-6] = 0
        vectors_low = vectors_high
        vectors_FmF_low = vectors_P32_FmF
        it += 1
    return energy_vec, vectors, vectors_P32_FmF, position_check_matrix

def generate_frequencies_at(b, f0 = 546.1206708):
    energy_lvl_D52, vectors, vectors_D52_FmF, p = Energy_at(b)
    energy_lvl_S12, vectors_S12, vectors_S12_FmF = S12_Energy_at(b)
    T_S12, Table_S12 = CG_mat_S12(0.5, 1.5)
    T_D52, Table_D52 = transformToFMFBasis(2.5, 1.5)
    
    energy_D52 = np.asarray(energy_lvl_D52).flatten()
    energy_S12 = np.asarray(np.diag(energy_lvl_S12)).flatten()
    # print(energy_D52, energy_S12)
    energys_table = np.zeros((24, 5))
    col_index = 4
    for i in range(8, 3, -1):
        row_index = 0
        for j in range(23, -1, -1):
            energys_table[row_index, col_index] = energy_D52[j] - energy_S12[i - 1]
            row_index += 1
        col_index -= 1
    
    h = 1
    freqs = np.abs(energys_table) / (h * 1e6)
    freqs_cali = freqs - (freqs[5, 2] - f0)
    
    num_rows = 24
    num_cols = 5
    Final_freqs = np.full((num_rows, num_cols), np.nan)
    
    col_labels = np.arange(-2, 3)
    row_labels = np.concatenate([
        np.flip(np.arange(-1, 2)),
        np.flip(np.arange(-2, 3)),
        np.flip(np.arange(-3, 4)),
        np.flip(np.arange(-4, 5))
    ])
    
    for row in range(num_rows):
        for col in range(num_cols):
            if abs(row_labels[row] - col_labels[col]) <= 2:
                Final_freqs[row, col] = freqs_cali[row, col]
    return Final_freqs

def fit_B_for_transition(f0, f1, row=22, col=1, B_min=0.0, B_max=10.0, tol=1e-9, max_iter=200):
    def transition_freq(B, f0=f0):
        freqs = generate_frequencies_at(B, f0)
        return freqs[row, col]

    def objective(B):
        return (transition_freq(B) - f1) ** 2

    res = minimize_scalar(
        objective,
        bounds=(B_min, B_max),
        method="bounded",
        options={"xatol": tol, "maxiter": max_iter},
    )

    B_best = res.x
    f_best = transition_freq(B_best)
    return B_best, f_best

# =============================================================================
# Transition Strength
# =============================================================================

def generateLabels_D52(I, J, T):
    dim_I = int(2 * I + 1)
    dim_J = int(2 * J + 1)
    dim_F = dim_I * dim_J
    col_names = []
    row_names = [''] * dim_F
    idx_FmF = 0
    for F in reversed(range(int(abs(J - I)), int(J + I) + 1)):
        for mF in reversed(range(-F, F + 1)):
            idx_FmF += 1
            col_names.append(f'{F},{mF}')
            idx_IJmImJ = 0
            for mJ in reversed(np.arange(-J, J + 1)):
                for mI in reversed(np.arange(-I, I + 1)):
                    idx_IJmImJ += 1
                    if idx_FmF == 1:
                        row_names[idx_IJmImJ - 1] = f'{mJ:.1f},{mI:.1f}'
    T_table = pd.DataFrame(T, index=row_names, columns=col_names)
    return col_names, row_names, T_table


def generateLabels_S12(I, J, T):
    dim_I = int(2 * I + 1)
    dim_J = int(2 * J + 1)
    dim_F = dim_I * dim_J
    col_names = []
    row_names = [''] * dim_F
    idx_FmF = 0
    F = 1
    for mF in range(-F, F + 1):
        idx_FmF += 1
        col_names.append(f'{F},{mF}')
    F = 2
    for mF in reversed(range(-F, F + 1)):
        idx_FmF += 1
        col_names.append(f'{F},{mF}')
        idx_IJmImJ = 0
        for mJ in reversed(np.arange(-J, J + 1)):
            for mI in reversed(np.arange(-I, I + 1)):
                idx_IJmImJ += 1
                if idx_FmF == 6:
                    row_names[idx_IJmImJ - 1] = f'{mJ:.1f},{mI:.1f}'
    T_table = pd.DataFrame(T, index=row_names, columns=col_names)
    return col_names, row_names, T_table


def generateLabels_P32(I, J, T):
    dim_I = int(2 * I + 1)
    dim_J = int(2 * J + 1)
    dim_F = dim_I * dim_J
    col_names = []
    row_names = [''] * dim_F
    idx_FmF = 0
    for F in reversed(range(int(abs(J - I)), int(J + I) + 1)):
        for mF in reversed(range(-F, F + 1)):
            idx_FmF += 1
            col_names.append(f'{F},{mF}')
            idx_IJmImJ = 0
            for mJ in reversed(np.arange(-J, J + 1)):
                for mI in reversed(np.arange(-I, I + 1)):
                    idx_IJmImJ += 1
                    if idx_FmF == 1:
                        row_names[idx_IJmImJ - 1] = f'{mJ:.1f},{mI:.1f}'
    T_table = pd.DataFrame(T, index=row_names, columns=col_names)
    return col_names, row_names, T_table

def parseLabel(label):
    a, b = label.split(',')
    return float(a), float(b)

def TransitionStrength(b, theta_k = 45, theta_p = 58):
    _, v_d52, _, _ = Energy_at(b)
    _, v_s12, _ = S12_Energy_at(b)
    tk = theta_k * np.pi / 180
    tp = theta_p * np.pi / 180
    G = np.array([
        0.25 * abs(np.cos(tp) * np.sin(2 * tk) - 2j * np.sin(tp) * np.sin(tk)),
        0.5  * abs(1j * np.sin(tp) * np.cos(tk) + np.cos(tp) * np.cos(2 * tk)),
        np.sqrt(6) / 4 * abs(np.cos(tp) * np.sin(2 * tk)),
        0.5  * abs(1j * np.sin(tp) * np.cos(tk) - np.cos(tp) * np.cos(2 * tk)),
        0.25 * abs(np.cos(tp) * np.sin(2 * tk) - 2j * np.sin(tp) * np.sin(tk))
    ])
    col_D52, row_D52, c_D52 = generateLabels_D52(2.5, 1.5, v_d52)
    col_S12, row_S12, c_S12 = generateLabels_S12(0.5, 1.5, v_s12)
    row_map_D = {s: i for i, s in enumerate(row_D52)}
    # print(c_S12)
    row_map_S = {s: i for i, s in enumerate(row_S12)}
    ts = np.zeros((len(col_D52), len(col_S12)))
    for i_d, lbl_d in enumerate(col_D52):
        F_D, mF_D = parseLabel(lbl_d)
        for i_s, lbl_s in enumerate(col_S12):
            F_S, mF_S = parseLabel(lbl_s)
            if abs(mF_D - mF_S) <= 2:
                s = 0
                for mJ_D in [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]:
                    for mI_D in [-1.5, -0.5, 0.5, 1.5]:
                        for mJ_S in [-0.5, 0.5]:
                            mI_S = mI_D
                            idx_D = row_map_D[f'{mI_D:.1f},{mJ_D:.1f}']
                            idx_S = row_map_S[f'{mI_S:.1f},{mJ_S:.1f}']
                            c_d = c_D52.iloc[idx_D, i_d]
                            c_s = c_S12.iloc[idx_S, i_s]
                            # print(idx_D, idx_S, c_d, c_s)
                            if c_d and c_s:
                                s += c_d * c_s * clebschgordan1(0.5, mJ_S, 2, mJ_D - mJ_S, 2.5, mJ_D)
                if s != 0:
                    ts[abs(i_d - 23), i_s] = abs(G[int(mF_D - mF_S + 2)] * s)
    return ts

def TransitionStrength_S12_P32(b, theta_k = 0, theta_p = 0):
    _, v_P32, _, _ = P32_Energy_at(b)
    _, v_s12, _ = S12_Energy_at(b)
    tk = theta_k * np.pi / 180
    tp = theta_p * np.pi / 180
    G = np.array([
        0 * abs(np.cos(tp) * np.sin(2 * tk) - 2j * np.sin(tp) * np.sin(tk)),
        1  * abs(1j * np.sin(tp) * np.cos(tk) + np.cos(tp) * np.cos(2 * tk)),
        np.sqrt(6) / 4 * abs(np.cos(tp) * np.sin(2 * tk)),
        1  * abs(1j * np.sin(tp) * np.cos(tk) - np.cos(tp) * np.cos(2 * tk)),
        0 * abs(np.cos(tp) * np.sin(2 * tk) - 2j * np.sin(tp) * np.sin(tk))
    ])
    G = np.array([0,1,1,1,0])
    col_P32, row_P32, c_P32 = generateLabels_P32(1.5, 1.5, v_P32)
    col_S12, row_S12, c_S12 = generateLabels_S12(0.5, 1.5, v_s12)
    row_map_P = {s: i for i, s in enumerate(row_P32)}
    # print(c_S12)
    row_map_S = {s: i for i, s in enumerate(row_S12)}
    ts = np.zeros((len(col_P32), len(col_S12)))
    for i_p, lbl_p in enumerate(col_P32):
        F_P, mF_P = parseLabel(lbl_p)
        for i_s, lbl_s in enumerate(col_S12):
            F_S, mF_S = parseLabel(lbl_s)
            if abs(mF_P - mF_S) <= 1:
                s = 0
                for mJ_P in [-1.5, -0.5, 0.5, 1.5]:
                    for mI_P in [-1.5, -0.5, 0.5, 1.5]:
                        for mJ_S in [-0.5, 0.5]:
                            mI_S = mI_P
                            idx_P = row_map_P[f'{mI_P:.1f},{mJ_P:.1f}']
                            idx_S = row_map_S[f'{mI_S:.1f},{mJ_S:.1f}']
                            c_p = c_P32.iloc[idx_P, i_p]
                            c_s = c_S12.iloc[idx_S, i_s]
                            # print(idx_P, idx_S, c_p, c_s)
                            if c_p and c_s:
                                s += c_p * c_s * clebschgordan1(0.5, mJ_S, 1, mJ_P - mJ_S, 1.5, mJ_P)
                if s != 0:
                    ts[abs(i_p - 15), i_s] = abs(G[int(mF_P - mF_S + 2)] * s)
    return ts

def TransitionStrength_D52_P32(b, theta_k = 0, theta_p = 0):
    _, v_d52, _, _ = Energy_at(b)
    _, v_p32, _, _ = P32_Energy_at(b)
    tk = theta_k * np.pi / 180
    tp = theta_p * np.pi / 180
    G = np.array([
        0.25 * abs(np.cos(tp) * np.sin(2 * tk) - 2j * np.sin(tp) * np.sin(tk)),
        0.5  * abs(1j * np.sin(tp) * np.cos(tk) + np.cos(tp) * np.cos(2 * tk)),
        np.sqrt(6) / 4 * abs(np.cos(tp) * np.sin(2 * tk)),
        0.5  * abs(1j * np.sin(tp) * np.cos(tk) - np.cos(tp) * np.cos(2 * tk)),
        0.25 * abs(np.cos(tp) * np.sin(2 * tk) - 2j * np.sin(tp) * np.sin(tk))
    ])
    G = np.array([0,1,0,1,0])
    col_D52, row_D52, c_D52 = generateLabels_D52(2.5, 1.5, v_d52)
    col_P32, row_P32, c_P32 = generateLabels_P32(1.5, 1.5, v_p32)
    row_map_D = {s: i for i, s in enumerate(row_D52)}
    # print(c_P32)
    row_map_P = {s: i for i, s in enumerate(row_P32)}
    ts = np.zeros((len(col_D52), len(col_P32)))
    for i_d, lbl_d in enumerate(col_D52):
        F_D, mF_D = parseLabel(lbl_d)
        for i_p, lbl_p in enumerate(col_P32):
            F_P, mF_P = parseLabel(lbl_p)
            if abs(mF_D - mF_P) <= 1:
                s = 0
                for mJ_D in [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]:
                    for mI_D in [-1.5, -0.5, 0.5, 1.5]:
                        for mJ_P in [-1.5, -0.5, 0.5, 1.5]:
                            mI_P = mI_D
                            idx_D = row_map_D[f'{mI_D:.1f},{mJ_D:.1f}']
                            idx_P = row_map_P[f'{mI_P:.1f},{mJ_P:.1f}']
                            c_d = c_D52.iloc[idx_D, i_d]
                            c_p = c_P32.iloc[idx_P, i_p]
                            # print(idx_D, idx_P, c_d, c_p)
                            if c_d and c_p:
                                s += c_d * c_p * clebschgordan1(1.5, mJ_P, 1, mJ_D - mJ_P, 2.5, mJ_D)
                if s != 0:
                    ts[abs(i_d - 23), i_p] = abs(G[int(mF_D - mF_P + 2)] * s)
    return ts

def TransitionStrength_D52_D52(b):
    _, v_d52, _, _ = Energy_at(b)
    col_D52, row_D52, c_D52 = generateLabels_D52(2.5, 1.5, v_d52)
    row_map_D = {s: i for i, s in enumerate(row_D52)}
    n = len(col_D52)
    ts = np.zeros((n, n))
    for i_i, lbl_i in enumerate(col_D52):
        _, mF_i = parseLabel(lbl_i)
        for i_f, lbl_f in enumerate(col_D52):
            _, mF_f = parseLabel(lbl_f)
            q = mF_f - mF_i
            if abs(q) > 1:
                continue
            amp = 0.0
            for mJ_i in [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]:
                mJ_f = mJ_i + q
                if not -2.5 <= mJ_f <= 2.5:
                    continue
                for mI in [-1.5, -0.5, 0.5, 1.5]:
                    idx_i = row_map_D[f'{mI:.1f},{mJ_i:.1f}']
                    idx_f = row_map_D[f'{mI:.1f},{mJ_f:.1f}']
                    c_i = c_D52.iloc[idx_i, i_i]
                    c_f = c_D52.iloc[idx_f, i_f]
                    if c_i == 0 or c_f == 0:
                        continue
                    amp += c_i * c_f * clebschgordan1(2.5, mJ_i, 1, q, 2.5, mJ_f)
            if amp != 0:
                ts[i_i, i_f] = abs(amp)
    return ts

# =============================================================================
# Magnetic Field Sensitivity
# =============================================================================

def Calculate_raw_freqs(b):
    h = 1
    # energy_lvl_D52 is already a 24-element vector from Energy_at
    energy_lvl_D52, _, _, _ = Energy_at(b)
    # energy_lvl_S12 is a diagonal matrix -> take its diagonal to get the 8-vector
    energy_lvl_S12, _, _ = S12_Energy_at(b)
    energy_D52 = energy_lvl_D52
    energy_S12 = np.diag(energy_lvl_S12)

    energys_table = np.zeros((24, 5))
    col_index = 4                      # MATLAB col_index = 5 -> zero-based = 4
    for i in range(7, 2, -1):          # fliplr(4:8)  -> indices 7,6,5,4
        row_index = 0
        for j in range(23, -1, -1):    # fliplr(1:24) -> 23 ... 0
            energys_table[row_index, col_index] = energy_D52[j] - energy_S12[i]
            row_index += 1
        col_index -= 1

    raw_freqs = np.abs(energys_table) / (h * 1e6)
    return raw_freqs

def generate_magnetic_field_sensitivity(b = 4.2095):
    delta_b = 0.0001
    raw_freqs_1 = Calculate_raw_freqs(b)
    raw_freqs_2 = Calculate_raw_freqs(b + delta_b)
    Sensitivity_matrix = (raw_freqs_1 - raw_freqs_2) / delta_b
    S = Sensitivity_matrix
    return S
