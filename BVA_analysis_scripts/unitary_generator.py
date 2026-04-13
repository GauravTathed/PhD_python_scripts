import numpy as np


def virtual_z_diagonal(z_frame):
    """Return the diagonal unitary corresponding to a virtual-Z frame."""
    return np.diag(np.exp(1j * np.asarray(z_frame, dtype=float)))


def ion_pulse_unitary(coupling, theta, phi, dim):
    """Return the embedded two-level rotation applied by one pulse."""
    i, j = (int(coupling[0]), int(coupling[1]))
    if not (0 <= i < dim and 0 <= j < dim):
        raise ValueError(f"coupling {coupling} is out of range for dim={dim}")
    if i == j:
        raise ValueError("coupling indices must be distinct")

    U = np.eye(dim, dtype=complex)
    c = np.cos(float(theta) / 2.0)
    s = np.sin(float(theta) / 2.0)
    phase = np.exp(1j * float(phi))

    U[i, i] = c
    U[j, j] = c
    U[i, j] = -1j * phase * s
    U[j, i] = -1j * np.conjugate(phase) * s
    return U


def infer_dimension(couplings, dim=None):
    """Infer the Hilbert-space dimension from the largest addressed level."""
    if dim is not None:
        dim = int(dim)
        if dim < 1:
            raise ValueError("dim must be at least 1")
        return dim

    normalized = [tuple(int(value) for value in coupling) for coupling in couplings]
    if len(normalized) == 0:
        raise ValueError("dim must be provided when couplings is empty")
    return max(max(coupling) for coupling in normalized) + 1


def unitary_from_pulse_sequence(couplings, thetas, phases, dim=None, final_frame=None, apply_final_frame=False):
    """
    Build the unitary implemented by a pulse sequence.

    Parameters
    ----------
    couplings : sequence of (i, j)
        Addressed two-level subspaces in application order.
    thetas : sequence of float
        Pulse angles in radians.
    phases : sequence of float
        Pulse phases in radians.
    dim : int, optional
        Hilbert-space dimension. If omitted, it is inferred from ``couplings``.
    final_frame : sequence of float, optional
        Per-level virtual-Z frame in radians.
    apply_final_frame : bool, optional
        If True, left-multiply the pulse sequence by ``diag(exp(1j * final_frame))``.

    Returns
    -------
    np.ndarray
        The implemented ``dim x dim`` unitary matrix.
    """
    normalized_couplings = [tuple(int(value) for value in coupling) for coupling in couplings]
    thetas = np.asarray(thetas, dtype=float)
    phases = np.asarray(phases, dtype=float)

    if len(normalized_couplings) != len(thetas) or len(normalized_couplings) != len(phases):
        raise ValueError("couplings, thetas, and phases must all have the same length")

    dim = infer_dimension(normalized_couplings, dim=dim)
    U = np.eye(dim, dtype=complex)
    for coupling, theta, phi in zip(normalized_couplings, thetas, phases):
        U = ion_pulse_unitary(coupling, theta, phi, dim) @ U

    if apply_final_frame:
        if final_frame is None:
            raise ValueError("final_frame must be provided when apply_final_frame=True")
        final_frame = np.asarray(final_frame, dtype=float)
        if final_frame.shape != (dim,):
            raise ValueError(f"final_frame must have shape ({dim},)")
        U = virtual_z_diagonal(final_frame) @ U

    return U


__all__ = [
    "virtual_z_diagonal",
    "ion_pulse_unitary",
    "infer_dimension",
    "unitary_from_pulse_sequence",
]
