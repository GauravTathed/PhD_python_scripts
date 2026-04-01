import numpy as np
from numpy.linalg import norm

from U_decomp import unitary_to_G_rotations


def wrap_pi(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def wrap_frame(z_frame):
    z = np.asarray(z_frame, dtype=float)
    return np.array([wrap_pi(value) for value in z], dtype=float)


def diagonal_phase_vector(V, tol=1e-8):
    V = np.asarray(V, dtype=complex)
    off_diag = V - np.diag(np.diag(V))
    if norm(off_diag) > tol:
        raise ValueError(f"Residual V is not diagonal enough: {norm(off_diag)}")

    diag = np.diag(V)
    if np.max(np.abs(np.abs(diag) - 1.0)) > 1e-7:
        raise ValueError("Residual V is diagonal but not unitary on the diagonal")
    return wrap_frame(np.angle(diag))


def rotation_matrix_to_pulse_parameters(G, tol=1e-10):
    G = np.asarray(G, dtype=complex)
    if G.ndim != 2 or G.shape[0] != G.shape[1]:
        raise ValueError("G must be square")

    dim = G.shape[0]
    candidates = [
        (p, q)
        for p in range(dim)
        for q in range(p + 1, dim)
        if abs(G[p, q]) > tol or abs(G[q, p]) > tol
    ]
    if len(candidates) != 1:
        raise ValueError(f"Expected exactly one coupled pair, found {candidates}")

    i, j = candidates[0]
    c = G[i, i]
    s = G[i, j]
    gamma = float(wrap_pi(np.angle(c)))
    theta = float(2.0 * np.arctan2(np.clip(np.abs(s), 0.0, 1.0), np.clip(np.abs(c), 0.0, 1.0)))
    phi = float((np.angle(s) - gamma + np.pi / 2) % (2 * np.pi))
    return {
        "coupling": (i, j),
        "theta": theta,
        "fraction": float(theta / np.pi),
        "phi": phi,
        "gamma": gamma,
    }


def minimize_gamma_for_step(step):
    alternate = dict(step)
    alternate["theta"] = float(2 * np.pi - step["theta"])
    alternate["fraction"] = float(alternate["theta"] / np.pi)
    alternate["phi"] = float((step["phi"] + np.pi) % (2 * np.pi))
    alternate["gamma"] = float(wrap_pi(step["gamma"] + np.pi))
    return alternate if abs(alternate["gamma"]) < abs(step["gamma"]) else dict(step)


def decompose_unitary_sequence(unitaries, center=0, tol=1e-10, minimize_gamma=True, initial_frame=None, return_metadata=False):
    if not unitaries:
        raise ValueError("Provide at least one unitary")

    unitary_list = [np.asarray(U, dtype=complex) for U in unitaries]
    dim = unitary_list[0].shape[0]
    for idx, U in enumerate(unitary_list, start=1):
        if U.ndim != 2 or U.shape[0] != U.shape[1]:
            raise ValueError(f"Unitary {idx} is not square")
        if U.shape[0] != dim:
            raise ValueError("All unitaries must have the same dimension")

    frame = np.zeros(dim, dtype=float) if initial_frame is None else np.asarray(initial_frame, dtype=float).copy()
    if frame.shape != (dim,):
        raise ValueError("initial_frame must have shape (dim,)")

    flat_schedule = []
    block_summaries = []

    for block_idx, U in enumerate(unitary_list, start=1):
        rotation_mats, V = unitary_to_G_rotations(U, center=center, tol=tol)
        V_phase = diagonal_phase_vector(V)

        incoming_frame = frame.copy()
        frame = wrap_frame(frame + V_phase)
        compile_frame_start = frame.copy()

        target_steps = [
            rotation_matrix_to_pulse_parameters(G.conj().T, tol=tol)
            for G in rotation_mats[::-1]
        ]
        if minimize_gamma:
            target_steps = [minimize_gamma_for_step(step) for step in target_steps]

        for step_idx, step in enumerate(target_steps, start=1):
            i, j = step["coupling"]
            phi_programmed = float((step["phi"] - (frame[i] - frame[j])) % (2 * np.pi))
            flat_schedule.append({
                "block": block_idx,
                "step_in_block": step_idx,
                "coupling": step["coupling"],
                "theta": step["theta"],
                "fraction": step["fraction"],
                "phi_programmed": phi_programmed,
                "phi_physical": step["phi"],
                "gamma": step["gamma"],
            })
            frame[i] += step["gamma"]
            frame[j] -= step["gamma"]
            frame = wrap_frame(frame)

        block_summaries.append({
            "block": block_idx,
            "incoming_frame": incoming_frame,
            "residual_from_V": V_phase,
            "compile_frame_start": compile_frame_start,
            "outgoing_frame": frame.copy(),
            "pulse_count": len(target_steps),
        })

    if return_metadata:
        return flat_schedule, frame.copy(), block_summaries
    return flat_schedule


def parse_string_to_int_list(s):
    s = str(s).strip()
    if s.startswith("[") and s.endswith("]"):
        inside = s[1:-1].strip()
        if not inside:
            return []
        return [int(value.strip()) for value in inside.split(",")]
    return [int(s)]


def build_transitions_list(index_to_str_map, couplings):
    missing = sorted({index for coupling in couplings for index in coupling if index not in index_to_str_map})
    if missing:
        raise KeyError(f"Mapping is missing indices: {missing}")

    transitions_list = []
    for start_idx, end_idx in couplings:
        left_list = parse_string_to_int_list(index_to_str_map[start_idx])
        right_list = parse_string_to_int_list(index_to_str_map[end_idx])
        transitions_list.append(left_list + right_list)
    return transitions_list


def qudit_qft(d, inverse=False):
    if d < 2:
        raise ValueError("d must be at least 2")
    omega = np.exp((-2j if inverse else 2j) * np.pi / d)
    j = np.arange(d).reshape((d, 1))
    k = np.arange(d).reshape((1, d))
    return omega ** (j * k) / np.sqrt(d)


def single_qudit_phase_oracle(hidden_value, d):
    hidden_value = int(hidden_value)
    if d < 2:
        raise ValueError("d must be at least 2")
    if not (0 <= hidden_value < d):
        raise ValueError("hidden_value must satisfy 0 <= hidden_value < dimension")
    omega = np.exp(2j * np.pi / d)
    phases = np.array([omega ** (hidden_value * x) for x in range(d)], dtype=complex)
    return np.diag(phases)


def grover_oracle(marked_state, d):
    marked_state = int(marked_state)
    if not (0 <= marked_state < d):
        raise ValueError("marked_state must satisfy 0 <= marked_state < dimension")
    phases = np.ones(d, dtype=complex)
    phases[marked_state] = -1.0
    return np.diag(phases)


def grover_diffuser(d):
    uniform_state = np.ones(d, dtype=complex) / np.sqrt(d)
    return 2.0 * np.outer(uniform_state, np.conjugate(uniform_state)) - np.eye(d, dtype=complex)


def build_algorithm_unitaries(algorithm_name, dimension, **algorithm_kwargs):
    algorithm_key = str(algorithm_name).strip().lower()
    dimension = int(dimension)
    if dimension < 2:
        raise ValueError("dimension must be at least 2")

    if algorithm_key in {"bva", "bv", "bernstein-vazirani", "bernstein_vazirani"}:
        if "hidden_value" not in algorithm_kwargs:
            raise ValueError("BVA requires hidden_value")
        hidden_value = algorithm_kwargs["hidden_value"]
        return [
            qudit_qft(dimension),
            single_qudit_phase_oracle(hidden_value, dimension),
            qudit_qft(dimension, inverse=True),
        ]

    if algorithm_key in {"grover", "grovers", "grover_search"}:
        if "marked_state" not in algorithm_kwargs:
            raise ValueError("Grover requires marked_state")
        marked_state = algorithm_kwargs["marked_state"]
        grover_iterations = int(algorithm_kwargs.get("grover_iterations", 1))
        if grover_iterations < 1:
            raise ValueError("grover_iterations must be at least 1")

        oracle = grover_oracle(marked_state, dimension)
        diffuser = grover_diffuser(dimension)
        return [qudit_qft(dimension)] + [gate for _ in range(grover_iterations) for gate in (oracle, diffuser)]

    raise ValueError(f"Unsupported algorithm_name: {algorithm_name}")


def compile_algorithm_pulses(mappings, dimension, algorithm_name, center=0, tol=1e-10, minimize_gamma=False, initial_frame=None, return_metadata=False, **algorithm_kwargs):
    """
    Compile a named algorithm into star-topology pulse lists.

    Parameters
    ----------
    mappings : dict[int, str]
        Map from logical level index to the transition label string used by the experiment.
        Each value should look like '0' or '[F, mF]' so it can be expanded into the pulse train format.
    dimension : int
        Qudit dimension.
    algorithm_name : str
        Currently supports 'BVA' and 'Grover'.
    return_metadata : bool, optional
        When True, also return a metadata dict with logical couplings, compiled schedule,
        final frame, and the source unitary list.

    Additional keyword arguments
    ----------------------------
    BVA:
        hidden_value : int
    Grover:s.py
        marked_state : int
        grover_iterations : int, optional

    Returns
    -------
    pulse_train_list, pulse_angle_list, pulse_phases_list
        By default.
    pulse_train_list, pulse_angle_list, pulse_phases_list, metadata
        When return_metadata=True.
    """
    unitaries = build_algorithm_unitaries(algorithm_name, dimension, **algorithm_kwargs)
    flat_schedule, final_frame, block_summaries = decompose_unitary_sequence(
        unitaries,
        center=center,
        tol=tol,
        minimize_gamma=minimize_gamma,
        initial_frame=initial_frame,
        return_metadata=True,
    )
    couplings = [step["coupling"] for step in flat_schedule]
    pulse_train_list = build_transitions_list(mappings, couplings)
    pulse_angle_list = [step["theta"] for step in flat_schedule]
    pulse_phases_list = [step["phi_programmed"] for step in flat_schedule]

    if not return_metadata:
        return pulse_train_list, pulse_angle_list, pulse_phases_list

    metadata = {
        "unitaries": unitaries,
        "flat_schedule": flat_schedule,
        "block_summaries": block_summaries,
        "couplings": couplings,
        "fractions": [step["fraction"] for step in flat_schedule],
        "gammas": [step["gamma"] for step in flat_schedule],
        "blocks": [step["block"] for step in flat_schedule],
        "steps_in_block": [step["step_in_block"] for step in flat_schedule],
        "final_frame": final_frame,
    }
    return pulse_train_list, pulse_angle_list, pulse_phases_list, metadata


__all__ = [
    "compile_algorithm_pulses",
    "build_algorithm_unitaries",
    "qudit_qft",
    "single_qudit_phase_oracle",
    "grover_oracle",
    "grover_diffuser",
]
