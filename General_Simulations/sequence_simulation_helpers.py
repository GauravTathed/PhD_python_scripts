import io
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.special import voigt_profile

try:
    import pandas as pd
except ImportError:
    pd = None

try:
    from IPython.display import display
except ImportError:
    display = None


ROW_LABELS = [[F, F - step] for F in [1, 2, 3, 4] for step in range(2 * F + 1)]
COL_LABELS = [-2, -1, 0, 1, 2]


def wrap_pi(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def wrap_frame(z_frame):
    z = np.asarray(z_frame, dtype=float)
    return np.array([wrap_pi(value) for value in z], dtype=float)


def virtual_z_diagonal(z_frame):
    return np.diag(np.exp(1j * np.asarray(z_frame, dtype=float)))


def state_populations(state):
    state = np.asarray(state, dtype=complex)
    return np.abs(state) ** 2


def basis_state(index, dimension):
    if not (0 <= int(index) < int(dimension)):
        raise ValueError("index must satisfy 0 <= index < dimension")
    state = np.zeros(int(dimension), dtype=complex)
    state[int(index)] = 1.0
    return state


def normalize_state(state):
    vector = np.asarray(state, dtype=complex)
    if vector.ndim != 1:
        raise ValueError("state must be one-dimensional")
    magnitude = np.linalg.norm(vector)
    if magnitude == 0.0:
        raise ValueError("state must have non-zero norm")
    return vector / magnitude


def coupling_key(coupling):
    if coupling is None or len(coupling) != 2:
        raise ValueError(f"Invalid coupling: {coupling}")
    return tuple(int(value) for value in coupling)


def parse_string_to_int_list(label):
    label = str(label).strip()
    if label.startswith("[") and label.endswith("]"):
        inside = label[1:-1].strip()
        if not inside:
            return []
        return [int(value.strip()) for value in inside.split(",")]
    return [int(label)]


def build_transitions_list(index_to_str_map, couplings):
    missing = sorted({index for coupling in couplings for index in coupling if index not in index_to_str_map})
    if missing:
        raise KeyError(f"mappings is missing indices: {missing}")
    transitions = []
    for start_idx, end_idx in couplings:
        left = parse_string_to_int_list(index_to_str_map[start_idx])
        right = parse_string_to_int_list(index_to_str_map[end_idx])
        transitions.append(left + right)
    return transitions


def format_state_label(index, mappings):
    return str(mappings.get(index, index))


def format_transition_label(transition):
    if len(transition) == 3:
        return f"{transition[0]} -> [{transition[1]}, {transition[2]}]"
    if len(transition) == 4:
        return f"[{transition[0]}, {transition[1]}] -> [{transition[2]}, {transition[3]}]"
    return str(list(transition))


def format_coupling_label(coupling, mappings):
    i, j = coupling_key(coupling)
    return f"{format_state_label(i, mappings)} <-> {format_state_label(j, mappings)}"


def extract_sequence_fields(sequence_steps):
    if not sequence_steps:
        raise ValueError("sequence_steps must contain at least one step")

    wait_duration_keys = (
        "duration_us",
        "wait_us",
        "wait_duration_us",
        "delay_us",
        "time_us",
        "wait_time_us",
        "free_evolution_us",
        "tau_us",
        "duration",
    )

    def resolve_step_kind(step):
        raw_kind = step.get("kind")
        if raw_kind is not None:
            kind = str(raw_kind).strip().lower()
            if kind in {"pulse", "drive"}:
                return "pulse"
            if kind in {"wait", "delay", "free", "free_evolution", "idle"}:
                return "wait"
            raise ValueError(f"Unsupported sequence step kind: {raw_kind}")
        if step.get("theta") is None and any(key in step for key in wait_duration_keys):
            return "wait"
        return "pulse"

    def resolve_wait_duration_us(step, step_index):
        duration_value = None
        for key in wait_duration_keys:
            if key in step:
                duration_value = step[key]
                break
        if duration_value is None:
            raise KeyError(
                f"sequence_steps[{step_index}] wait step is missing a duration key. "
                f"Supported keys: {wait_duration_keys}"
            )
        duration_us = float(duration_value)
        if duration_us < 0.0:
            raise ValueError(f"sequence_steps[{step_index}] wait duration must be non-negative")
        return duration_us

    normalized_steps = []
    for pulse_index, step in enumerate(sequence_steps):
        if not isinstance(step, dict):
            raise TypeError(f"sequence_steps[{pulse_index}] must be a dict")
        kind = resolve_step_kind(step)
        coupling = step.get("coupling", step.get("states"))
        if coupling is None:
            raise KeyError(f"sequence_steps[{pulse_index}] is missing 'coupling'")
        normalized = {
            "kind": kind,
            "coupling": coupling_key(coupling),
            "theta": 0.0,
            "programmed_phase": 0.0,
            "duration_us": None,
            "pi_time_us": None if step.get("pi_time_us") is None else float(step.get("pi_time_us")),
            "sensitivity_mhz_per_g": None if step.get("sensitivity_mhz_per_g") is None else float(step.get("sensitivity_mhz_per_g")),
        }
        if kind == "pulse":
            theta = step.get("theta")
            phase = step.get("phi_programmed", step.get("phi", step.get("phase")))
            if theta is None:
                raise KeyError(f"sequence_steps[{pulse_index}] pulse step is missing 'theta'")
            if phase is None:
                raise KeyError(f"sequence_steps[{pulse_index}] pulse step is missing 'phi_programmed' or 'phi'")
            normalized["theta"] = float(theta)
            normalized["programmed_phase"] = float(phase)
        else:
            normalized["duration_us"] = resolve_wait_duration_us(step, pulse_index)
        normalized_steps.append(normalized)
    return normalized_steps


def display_records(records):
    if not records:
        print("No rows to display.")
        return
    if pd is not None and display is not None:
        display(pd.DataFrame(records))
        return
    headers = list(records[0].keys())
    print(" | ".join(headers))
    for record in records:
        print(" | ".join(str(record[key]) for key in headers))


@lru_cache(maxsize=None)
def load_csv_matrix(filename):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    matrix = np.loadtxt(path, delimiter=",")
    matrix = np.asarray(matrix, dtype=float)
    matrix[matrix == 0] = np.nan
    return matrix


def compute_pi_time_matrix(transition_strengths, reference_pi_times_us, reference_strength_indices):
    transition_strengths = np.asarray(transition_strengths, dtype=float)
    reference_pi_times_us = np.asarray(reference_pi_times_us, dtype=float)
    reference_strengths = np.array(
        [transition_strengths[row, col] for row, col in reference_strength_indices],
        dtype=float,
    )
    if not np.all(np.isfinite(reference_strengths)):
        raise ValueError("Reference strengths contain invalid entries")

    factors = reference_pi_times_us * reference_strengths
    pi_times = np.full_like(transition_strengths, np.nan, dtype=float)
    for row_index, row_label in enumerate(ROW_LABELS):
        for col_index in range(transition_strengths.shape[1]):
            strength = transition_strengths[row_index, col_index]
            if not np.isfinite(strength):
                continue
            delta_m = (row_label[1] - COL_LABELS[col_index]) + 2
            if 0 <= delta_m < len(factors):
                pi_times[row_index, col_index] = factors[delta_m] / strength
    return pi_times


def values_from_matrix(transitions, matrix):
    matrix = np.asarray(matrix, dtype=float)
    values = []
    for transition in transitions:
        row_label = [transition[1], transition[2]]
        col_label = transition[0]
        row_index = next((idx for idx, label in enumerate(ROW_LABELS) if label == row_label), None)
        col_index = COL_LABELS.index(col_label) if col_label in COL_LABELS else None
        if row_index is None or col_index is None:
            raise ValueError(f"Could not resolve transition {transition} in the calibration matrix")
        value = matrix[row_index, col_index]
        if not np.isfinite(value):
            raise ValueError(f"Matrix lookup returned an invalid value for transition {transition}")
        values.append(float(value))
    return np.asarray(values, dtype=float)


def resolve_pi_times_us(couplings, default_pi_time_us, pi_time_by_coupling_us=None, per_pulse_pi_time_us=None):
    couplings = [coupling_key(coupling) for coupling in couplings]
    if per_pulse_pi_time_us is not None:
        pi_times_us = np.asarray(per_pulse_pi_time_us, dtype=float)
        if pi_times_us.shape != (len(couplings),):
            raise ValueError("per_pulse_pi_time_us must have length equal to the pulse count")
        return pi_times_us
    default_pi_time_us = float(default_pi_time_us)
    pi_time_by_coupling_us = {} if pi_time_by_coupling_us is None else {
        coupling_key(key): float(value) for key, value in pi_time_by_coupling_us.items()
    }
    return np.asarray([pi_time_by_coupling_us.get(coupling, default_pi_time_us) for coupling in couplings], dtype=float)


def pulse_durations_us_from_thetas(thetas, pi_times_us):
    thetas = np.asarray(thetas, dtype=float)
    pi_times_us = np.asarray(pi_times_us, dtype=float)
    if thetas.shape != pi_times_us.shape:
        raise ValueError("thetas and pi_times_us must have the same shape")
    return thetas * pi_times_us / np.pi


def pulse_timing_us(durations_us):
    durations_us = np.asarray(durations_us, dtype=float)
    t_end = np.cumsum(durations_us)
    t_start = np.concatenate(([0.0], t_end[:-1])) if len(durations_us) else np.asarray([], dtype=float)
    return np.column_stack((t_start, t_end)) if len(durations_us) else np.empty((0, 2), dtype=float)


def load_line_fit_parameters(filename):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"Line-signal fit file not found: {path}")

    B0_fit = None
    B0_err = None
    harmonic_lines = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if B0_fit is None:
                values = np.fromstring(stripped, sep=" ")
                B0_fit, B0_err = values[0], values[1]
            else:
                harmonic_lines.append(stripped)
    harmonics = np.loadtxt(io.StringIO("\n".join(harmonic_lines)))
    if harmonics.ndim == 1:
        harmonics = harmonics[None, :]
    freqs = harmonics[:, 0]
    amplitudes = harmonics[:, 1]
    phases = harmonics[:, 3]
    return B0_fit, B0_err, freqs, amplitudes, phases


def default_line_noise_model(param_file=None):
    return {
        "B0_mG": 0.2289,
        "A60_mG": -0.2791,
        "phi60_rad": 1.0002,
        "A180_mG": -0.0774,
        "phi180_rad": -0.2650,
        "param_file": None if param_file is None else str(param_file),
    }


def _get_line_params(param_file, B0, A60, phi60, A180, phi180):
    if param_file is not None and Path(param_file).exists():
        B0_fit, _, freqs, amplitudes, phases = load_line_fit_parameters(param_file)
        return float(B0_fit), np.asarray(freqs, dtype=float), np.asarray(amplitudes, dtype=float), np.asarray(phases, dtype=float)
    freqs = np.array([60.0, 180.0], dtype=float)
    amplitudes = np.array([A60, A180], dtype=float)
    phases = np.array([phi60, phi180], dtype=float)
    return float(B0), freqs, amplitudes, phases


def _line_model_arrays(line_model):
    if line_model is None:
        line_model = default_line_noise_model()
    return _get_line_params(
        line_model.get("param_file"),
        line_model.get("B0_mG", 0.2289),
        line_model.get("A60_mG", -0.2791),
        line_model.get("phi60_rad", 1.0002),
        line_model.get("A180_mG", -0.0774),
        line_model.get("phi180_rad", -0.2650),
    )


def line_field_relative_G(t_sec, line_model=None):
    t_sec = np.asarray(t_sec, dtype=float)
    B0_use, freqs_hz, amplitudes_mG, phases_rad = _line_model_arrays(line_model)
    field_mG = np.full_like(t_sec, B0_use, dtype=float)
    reference_mG = B0_use
    for freq_hz, amplitude_mG, phase_rad in zip(freqs_hz, amplitudes_mG, phases_rad):
        field_mG += amplitude_mG * np.cos(2.0 * np.pi * freq_hz * t_sec + phase_rad)
        reference_mG += amplitude_mG * np.cos(phase_rad)
    return (field_mG - reference_mG) * 1e-3


def line_field_primitive_Gs(t_sec, line_model=None):
    t_sec = np.asarray(t_sec, dtype=float)
    B0_use, freqs_hz, amplitudes_mG, phases_rad = _line_model_arrays(line_model)
    B_line_0 = B0_use + np.sum(amplitudes_mG * np.cos(phases_rad))
    primitive_mGs = (B0_use - B_line_0) * t_sec
    for freq_hz, amplitude_mG, phase_rad in zip(freqs_hz, amplitudes_mG, phases_rad):
        omega = 2.0 * np.pi * freq_hz
        primitive_mGs += (amplitude_mG / omega) * (np.sin(omega * t_sec + phase_rad) - np.sin(phase_rad))
    return primitive_mGs * 1e-3


def compute_line_phase_corrections(durations_us, sensitivities_mhz_per_g, line_model=None, t_ref="start", enabled=True):
    durations_us = np.asarray(durations_us, dtype=float)
    sensitivities_mhz_per_g = np.asarray(sensitivities_mhz_per_g, dtype=float)
    if durations_us.shape != sensitivities_mhz_per_g.shape:
        raise ValueError("durations_us and sensitivities_mhz_per_g must have the same shape")
    if not enabled or len(durations_us) == 0:
        zeros = np.zeros(len(durations_us), dtype=float)
        return zeros, zeros, zeros

    timing_us = pulse_timing_us(durations_us)
    pulse_starts_sec = timing_us[:, 0] * 1e-6
    pulse_ends_sec = timing_us[:, 1] * 1e-6
    if t_ref == "start":
        t_refs_sec = pulse_starts_sec
    elif t_ref == "end":
        t_refs_sec = pulse_ends_sec
    elif t_ref == "mid":
        t_refs_sec = 0.5 * (pulse_starts_sec + pulse_ends_sec)
    else:
        raise ValueError("t_ref must be one of 'start', 'mid', or 'end'")

    area_Gs = line_field_primitive_Gs(pulse_starts_sec, line_model=line_model)
    B_rel_G = line_field_relative_G(t_refs_sec, line_model=line_model)
    phases_rad = 2.0 * np.pi * 1e6 * sensitivities_mhz_per_g * area_Gs
    detunings_mhz = B_rel_G * sensitivities_mhz_per_g
    return phases_rad, detunings_mhz, t_refs_sec


def laser_noise_widths_hz(gaussian_T_s, lorentzian_L_s):
    sigma_hz = np.sqrt(2.0) / float(gaussian_T_s)
    gamma_hz = 1.0 / float(lorentzian_L_s)
    return sigma_hz, gamma_hz


@lru_cache(maxsize=None)
def _cached_voigt_cdf_grid(sigma_hz, gamma_hz, x_min_hz, x_max_hz, num_points):
    sigma_hz = float(sigma_hz)
    gamma_hz = float(gamma_hz)
    x_min_hz = float(x_min_hz)
    x_max_hz = float(x_max_hz)
    num_points = int(num_points)
    if num_points < 2:
        raise ValueError("num_points must be at least 2")
    if not x_min_hz < x_max_hz:
        raise ValueError("x_min_hz must be less than x_max_hz")

    x_vals = np.linspace(x_min_hz, x_max_hz, num_points, dtype=float)
    dx = x_vals[1] - x_vals[0]
    profile_vals = voigt_profile(x_vals, sigma_hz, gamma_hz)
    cdf = np.cumsum(profile_vals) * dx
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0.0:
        raise ValueError("Voigt CDF normalization failed")
    cdf /= cdf[-1]
    cdf = np.maximum.accumulate(cdf)
    cdf[-1] = 1.0
    return x_vals, cdf


def sample_voigt_hz(rng=None, sigma_hz=0.0, gamma_hz=0.0, size=None, x_min_hz=-1000.0, x_max_hz=1000.0, num_points=100000):
    generator = np.random.default_rng() if rng is None else rng
    sigma_hz = max(float(sigma_hz), 0.0)
    gamma_hz = max(float(gamma_hz), 0.0)
    if sigma_hz == 0.0 and gamma_hz == 0.0:
        if size is None:
            return 0.0
        return np.zeros(size, dtype=float)
    x_vals, cdf = _cached_voigt_cdf_grid(sigma_hz, gamma_hz, x_min_hz, x_max_hz, num_points)
    samples = np.interp(generator.random(size=size), cdf, x_vals)
    if size is None:
        return float(samples)
    return samples


def sample_frequency_calibration_error_hz(rng=None, sigma_hz=0.0, hwhm_hz=14.976, size=None, x_min_hz=-1000.0, x_max_hz=1000.0, num_points=100000):
    return sample_voigt_hz(
        rng=rng,
        sigma_hz=sigma_hz,
        gamma_hz=float(hwhm_hz),
        size=size,
        x_min_hz=x_min_hz,
        x_max_hz=x_max_hz,
        num_points=num_points,
    )


def sample_laser_noise_hz(rng=None, gaussian_T_s=2 * np.pi * 2206.8675860990334e-6, lorentzian_L_s=2 * np.pi * 3353.155420102816e-6, x_min_hz=-1000.0, x_max_hz=1000.0, num_points=100000):
    sigma_hz, gamma_hz = laser_noise_widths_hz(gaussian_T_s, lorentzian_L_s)
    return sample_voigt_hz(
        rng=rng,
        sigma_hz=sigma_hz,
        gamma_hz=gamma_hz,
        size=None,
        x_min_hz=x_min_hz,
        x_max_hz=x_max_hz,
        num_points=num_points,
    )


def sample_rabi_frequency_scale(rng=None, std_fraction=0.037):
    generator = np.random.default_rng() if rng is None else rng
    return 1.0 + float(std_fraction) * generator.normal()


def sample_line_signal_field_offset_g(rng=None, std_g=50e-6):
    generator = np.random.default_rng() if rng is None else rng
    return float(std_g) * generator.normal()


def line_signal_field_noise_corrections(durations_us, sensitivities_mhz_per_g, line_field_offset_g=0.0, t_ref="start"):
    durations_us = np.asarray(durations_us, dtype=float)
    sensitivities_mhz_per_g = np.asarray(sensitivities_mhz_per_g, dtype=float)
    if durations_us.shape != sensitivities_mhz_per_g.shape:
        raise ValueError("durations_us and sensitivities_mhz_per_g must have the same shape")
    if len(durations_us) == 0 or float(line_field_offset_g) == 0.0:
        zeros = np.zeros(len(durations_us), dtype=float)
        return zeros, zeros, zeros

    timing_us = pulse_timing_us(durations_us)
    pulse_starts_sec = timing_us[:, 0] * 1e-6
    pulse_ends_sec = timing_us[:, 1] * 1e-6
    if t_ref == "start":
        t_refs_sec = pulse_starts_sec
    elif t_ref == "end":
        t_refs_sec = pulse_ends_sec
    elif t_ref == "mid":
        t_refs_sec = 0.5 * (pulse_starts_sec + pulse_ends_sec)
    else:
        raise ValueError("t_ref must be one of 'start', 'mid', or 'end'")

    area_Gs = float(line_field_offset_g) * t_refs_sec
    phases_rad = 2.0 * np.pi * 1e6 * sensitivities_mhz_per_g * area_Gs
    detunings_mhz = float(line_field_offset_g) * sensitivities_mhz_per_g
    return phases_rad, detunings_mhz, t_refs_sec


def build_unique_coupling_noise_inputs(couplings, pulse_sensitivities_mhz_per_g, atol=1e-12):
    couplings = [coupling_key(coupling) for coupling in couplings]
    pulse_sensitivities_mhz_per_g = np.asarray(pulse_sensitivities_mhz_per_g, dtype=float)
    if pulse_sensitivities_mhz_per_g.shape != (len(couplings),):
        raise ValueError("pulse_sensitivities_mhz_per_g must have length equal to the pulse count")
    unique_couplings, unique_sensitivities, index_lookup = [], [], {}
    for coupling, sensitivity in zip(couplings, pulse_sensitivities_mhz_per_g):
        if coupling not in index_lookup:
            index_lookup[coupling] = len(unique_couplings)
            unique_couplings.append(coupling)
            unique_sensitivities.append(float(sensitivity))
        else:
            previous = unique_sensitivities[index_lookup[coupling]]
            if not np.isclose(previous, float(sensitivity), atol=atol, rtol=0.0):
                raise ValueError(
                    f"Coupling {coupling} appears with multiple sensitivities: {previous} and {float(sensitivity)}"
                )
    return {
        "couplings": unique_couplings,
        "sensitivities_mhz_per_g": np.asarray(unique_sensitivities, dtype=float),
    }


def sample_unique_coupling_detunings_mhz(unique_sensitivities_mhz_per_g, rng=None, calibration_scale=1.0, magnetic_scale=1.0, laser_scale=1.0, calibration_noise_hwhm_hz=14.976, calibration_noise_sigma_hz=0.0, magnetic_field_std_g=11.13e-6, gaussian_T_s=2 * np.pi * 2206.8675860990334e-6, lorentzian_L_s=2 * np.pi * 3353.155420102816e-6, voigt_x_min_hz=-1000.0, voigt_x_max_hz=1000.0, voigt_num_points=100000):
    generator = np.random.default_rng() if rng is None else rng
    unique_sensitivities_mhz_per_g = np.asarray(unique_sensitivities_mhz_per_g, dtype=float)
    num_couplings = len(unique_sensitivities_mhz_per_g)

    shared_calibration_noise_hz = calibration_scale * sample_frequency_calibration_error_hz(
        generator,
        sigma_hz=calibration_noise_sigma_hz,
        hwhm_hz=calibration_noise_hwhm_hz,
        size=None,
        x_min_hz=voigt_x_min_hz,
        x_max_hz=voigt_x_max_hz,
        num_points=voigt_num_points,
    )
    calibration_noise_mhz = np.full(num_couplings, shared_calibration_noise_hz * 1e-6, dtype=float)

    magnetic_field_sample_g = magnetic_scale * generator.normal(0.0, magnetic_field_std_g)
    magnetic_noise_mhz = magnetic_field_sample_g * unique_sensitivities_mhz_per_g

    shared_laser_noise_hz = laser_scale * sample_laser_noise_hz(
        generator,
        gaussian_T_s=gaussian_T_s,
        lorentzian_L_s=lorentzian_L_s,
        x_min_hz=voigt_x_min_hz,
        x_max_hz=voigt_x_max_hz,
        num_points=voigt_num_points,
    )
    laser_noise_mhz = np.full(num_couplings, shared_laser_noise_hz * 1e-6, dtype=float)
    return calibration_noise_mhz + magnetic_noise_mhz + laser_noise_mhz


def simulate_schedule_state(couplings, durations_us, thetas, programmed_phases, initial_state, final_frame=None, apply_final_frame=False, line_phase_offsets_rad=None, line_detunings_mhz=None, global_detuning_mhz_by_coupling=None, rabi_frequency_scales=None):
    couplings = [coupling_key(coupling) for coupling in couplings]
    durations_us = np.asarray(durations_us, dtype=float)
    thetas = np.asarray(thetas, dtype=float)
    programmed_phases = np.asarray(programmed_phases, dtype=float)
    initial_state = np.asarray(initial_state, dtype=complex)
    num_pulses = len(couplings)
    if durations_us.shape != (num_pulses,) or thetas.shape != (num_pulses,) or programmed_phases.shape != (num_pulses,):
        raise ValueError("couplings, durations_us, thetas, and programmed_phases must all have the same length")

    line_phase_offsets_rad = np.zeros(num_pulses, dtype=float) if line_phase_offsets_rad is None else np.asarray(line_phase_offsets_rad, dtype=float)
    line_detunings_mhz = np.zeros(num_pulses, dtype=float) if line_detunings_mhz is None else np.asarray(line_detunings_mhz, dtype=float)
    rabi_frequency_scales = np.ones(num_pulses, dtype=float) if rabi_frequency_scales is None else np.asarray(rabi_frequency_scales, dtype=float)
    if line_phase_offsets_rad.shape != (num_pulses,) or line_detunings_mhz.shape != (num_pulses,) or rabi_frequency_scales.shape != (num_pulses,):
        raise ValueError("line-phase, line-detuning, and rabi-scale arrays must have length equal to the pulse count")

    dimension = initial_state.shape[0]
    state = initial_state.copy()
    H_det_global = np.zeros((dimension, dimension), dtype=complex)
    if global_detuning_mhz_by_coupling is not None:
        for coupling, delta_mhz in global_detuning_mhz_by_coupling.items():
            i, j = coupling_key(coupling)
            delta_rad_per_us = 2.0 * np.pi * float(delta_mhz)
            H_det_global[i, i] += 0.5 * delta_rad_per_us
            H_det_global[j, j] -= 0.5 * delta_rad_per_us

    total_phases = programmed_phases + line_phase_offsets_rad
    for (i, j), duration_us, theta, total_phase, delta_line_mhz, rabi_scale in zip(
        couplings,
        durations_us,
        thetas,
        total_phases,
        line_detunings_mhz,
        rabi_frequency_scales,
    ):
        if duration_us <= 0.0:
            continue
        omega_rad_per_us = float(rabi_scale) * theta / duration_us
        H_pulse = H_det_global.copy()
        delta_line_rad_per_us = 2.0 * np.pi * float(delta_line_mhz)
        H_pulse[i, i] += 0.5 * delta_line_rad_per_us
        H_pulse[j, j] -= 0.5 * delta_line_rad_per_us
        phase_factor = np.exp(1j * total_phase)
        H_pulse[i, j] += 0.5 * omega_rad_per_us * phase_factor
        H_pulse[j, i] += 0.5 * omega_rad_per_us * np.conjugate(phase_factor)
        state = expm(-1j * H_pulse * duration_us) @ state

    if apply_final_frame and final_frame is not None:
        state = virtual_z_diagonal(final_frame) @ state
    return state


def deterministic_sequence_result(couplings, thetas, programmed_phases, initial_state, pulse_sensitivities_mhz_per_g, pi_times_us=None, durations_us=None, final_frame=None, apply_final_frame=False, use_line_signal=True, line_model=None, line_reference="start", target_state=None):
    if durations_us is None:
        if pi_times_us is None:
            raise ValueError("Provide either durations_us or pi_times_us")
        durations_us = pulse_durations_us_from_thetas(thetas, pi_times_us)
    durations_us = np.asarray(durations_us, dtype=float)
    timing_us = pulse_timing_us(durations_us)
    line_phase_offsets_rad, line_detunings_mhz, line_reference_times_s = compute_line_phase_corrections(
        durations_us,
        pulse_sensitivities_mhz_per_g,
        line_model=line_model,
        t_ref=line_reference,
        enabled=use_line_signal,
    )
    final_state = simulate_schedule_state(
        couplings,
        durations_us,
        thetas,
        programmed_phases,
        initial_state=initial_state,
        final_frame=final_frame,
        apply_final_frame=apply_final_frame,
        line_phase_offsets_rad=line_phase_offsets_rad,
        line_detunings_mhz=line_detunings_mhz,
        global_detuning_mhz_by_coupling={coupling_key(coupling): 0.0 for coupling in couplings},
        rabi_frequency_scales=np.ones(len(couplings), dtype=float),
    )
    populations = state_populations(final_state)
    fidelity = np.nan
    if target_state is not None:
        fidelity = float(np.real_if_close(abs(np.vdot(normalize_state(target_state), final_state)) ** 2))
    return {
        "pi_times_us": None if pi_times_us is None else np.asarray(pi_times_us, dtype=float),
        "pulse_durations_us": durations_us,
        "durations_us": durations_us,
        "timing_us": timing_us,
        "line_phase_offsets_rad": line_phase_offsets_rad,
        "line_detunings_mhz": line_detunings_mhz,
        "line_reference_times_s": line_reference_times_s,
        "final_state": final_state,
        "populations": populations,
        "fidelity": fidelity,
    }


def monte_carlo_sequence_result(couplings, thetas, programmed_phases, initial_state, pulse_sensitivities_mhz_per_g, pi_times_us=None, durations_us=None, final_frame=None, apply_final_frame=True, use_line_signal=False, line_model=None, line_reference="start", num_shots=100, rng_seed=None, calibration_scale=1.0, magnetic_scale=1.0, laser_scale=1.0, calibration_noise_hwhm_hz=14.976, calibration_noise_sigma_hz=0.0, magnetic_field_std_g=11.13e-6, laser_gaussian_T_s=2 * np.pi * 2206.8675860990334e-6, laser_lorentzian_L_s=2 * np.pi * 3353.155420102816e-6, voigt_x_min_hz=-1000.0, voigt_x_max_hz=1000.0, voigt_num_points=100000, rabi_frequency_std_fraction=0.037, line_signal_noise_std_g=50e-6, target_state=None):
    couplings = [coupling_key(coupling) for coupling in couplings]
    thetas = np.asarray(thetas, dtype=float)
    programmed_phases = np.asarray(programmed_phases, dtype=float)
    pulse_sensitivities_mhz_per_g = np.asarray(pulse_sensitivities_mhz_per_g, dtype=float)
    initial_state = np.asarray(initial_state, dtype=complex)
    num_shots = int(num_shots)
    if num_shots < 1:
        raise ValueError("num_shots must be at least 1")

    if durations_us is None:
        if pi_times_us is None:
            raise ValueError("Provide either durations_us or pi_times_us")
        pi_times_us = np.asarray(pi_times_us, dtype=float)
        durations_us = pulse_durations_us_from_thetas(thetas, pi_times_us)
    else:
        durations_us = np.asarray(durations_us, dtype=float)
        pi_times_us = None if pi_times_us is None else np.asarray(pi_times_us, dtype=float)
    timing_us = pulse_timing_us(durations_us)
    base_line_phase_offsets_rad, base_line_detunings_mhz, line_reference_times_s = compute_line_phase_corrections(
        durations_us,
        pulse_sensitivities_mhz_per_g,
        line_model=line_model,
        t_ref=line_reference,
        enabled=use_line_signal,
    )
    unique_noise_inputs = build_unique_coupling_noise_inputs(couplings, pulse_sensitivities_mhz_per_g)
    generator = np.random.default_rng(rng_seed)

    noiseless_state = simulate_schedule_state(
        couplings,
        durations_us,
        thetas,
        programmed_phases,
        initial_state=initial_state,
        final_frame=final_frame,
        apply_final_frame=apply_final_frame,
        line_phase_offsets_rad=base_line_phase_offsets_rad,
        line_detunings_mhz=base_line_detunings_mhz,
        global_detuning_mhz_by_coupling={coupling: 0.0 for coupling in unique_noise_inputs["couplings"]},
        rabi_frequency_scales=np.ones(len(couplings), dtype=float),
    )

    target_state_normalized = normalize_state(target_state) if target_state is not None else None
    shot_populations, shot_fidelities = [], []
    for _ in range(num_shots):
        sampled_noise = sample_unique_coupling_detunings_mhz(
            unique_noise_inputs["sensitivities_mhz_per_g"],
            rng=generator,
            calibration_scale=calibration_scale,
            magnetic_scale=magnetic_scale,
            laser_scale=laser_scale,
            calibration_noise_hwhm_hz=calibration_noise_hwhm_hz,
            calibration_noise_sigma_hz=calibration_noise_sigma_hz,
            magnetic_field_std_g=magnetic_field_std_g,
            gaussian_T_s=laser_gaussian_T_s,
            lorentzian_L_s=laser_lorentzian_L_s,
            voigt_x_min_hz=voigt_x_min_hz,
            voigt_x_max_hz=voigt_x_max_hz,
            voigt_num_points=voigt_num_points,
        )
        detuning_map = {coupling: detuning_mhz for coupling, detuning_mhz in zip(unique_noise_inputs["couplings"], sampled_noise)}
        shared_rabi_scale = sample_rabi_frequency_scale(generator, std_fraction=rabi_frequency_std_fraction)
        rabi_frequency_scales = np.full(len(couplings), shared_rabi_scale, dtype=float)
        line_field_offset_g = sample_line_signal_field_offset_g(generator, std_g=line_signal_noise_std_g)
        line_noise_phases_rad, line_noise_detunings_mhz, _ = line_signal_field_noise_corrections(
            durations_us,
            pulse_sensitivities_mhz_per_g,
            line_field_offset_g=line_field_offset_g,
            t_ref=line_reference,
        )

        final_state = simulate_schedule_state(
            couplings,
            durations_us,
            thetas,
            programmed_phases,
            initial_state=initial_state,
            final_frame=final_frame,
            apply_final_frame=apply_final_frame,
            line_phase_offsets_rad=base_line_phase_offsets_rad + line_noise_phases_rad,
            line_detunings_mhz=base_line_detunings_mhz + line_noise_detunings_mhz,
            global_detuning_mhz_by_coupling=detuning_map,
            rabi_frequency_scales=rabi_frequency_scales,
        )
        shot_populations.append(state_populations(final_state))
        if target_state_normalized is not None:
            shot_fidelities.append(float(np.real_if_close(abs(np.vdot(target_state_normalized, final_state)) ** 2)))

    shot_populations = np.asarray(shot_populations, dtype=float)
    mean_populations = np.mean(shot_populations, axis=0)
    std_populations = np.std(shot_populations, axis=0, ddof=1) if num_shots > 1 else np.zeros_like(mean_populations)
    shot_fidelities = np.asarray(shot_fidelities, dtype=float)
    mean_fidelity = float(np.mean(shot_fidelities)) if len(shot_fidelities) else np.nan
    std_fidelity = float(np.std(shot_fidelities, ddof=1)) if len(shot_fidelities) > 1 else (0.0 if len(shot_fidelities) == 1 else np.nan)
    return {
        "pi_times_us": pi_times_us,
        "pulse_durations_us": durations_us,
        "durations_us": durations_us,
        "timing_us": timing_us,
        "line_phase_offsets_rad": base_line_phase_offsets_rad,
        "line_detunings_mhz": base_line_detunings_mhz,
        "line_reference_times_s": line_reference_times_s,
        "noiseless_state": noiseless_state,
        "mean_populations": mean_populations,
        "std_populations": std_populations,
        "shot_populations": shot_populations,
        "mean_fidelity": mean_fidelity,
        "std_fidelity": std_fidelity,
    }


def prepare_sequence_data(sequence_steps, dimension, mappings, default_pi_time_us=20.0, pi_time_by_coupling_us=None, per_pulse_pi_time_us=None, pulse_sensitivities_mhz_per_g=None, transition_strengths_file=None, sensitivity_matrix_file=None, reference_pi_times_us=None, reference_strength_indices=None):
    dimension = int(dimension)
    normalized_steps = extract_sequence_fields(sequence_steps)
    couplings = [step["coupling"] for step in normalized_steps]
    thetas = np.asarray([step["theta"] for step in normalized_steps], dtype=float)
    phases = np.asarray([step["programmed_phase"] for step in normalized_steps], dtype=float)
    step_kinds = [step["kind"] for step in normalized_steps]
    pulse_indices = [index for index, step in enumerate(normalized_steps) if step["kind"] == "pulse"]
    wait_indices = [index for index, step in enumerate(normalized_steps) if step["kind"] == "wait"]
    max_index = max(max(coupling) for coupling in couplings)
    if max_index >= dimension:
        raise ValueError(f"dimension={dimension} is too small for the largest coupling index {max_index}")

    transitions = build_transitions_list(mappings, couplings)
    pulse_couplings = [couplings[index] for index in pulse_indices]
    pulse_transitions = [transitions[index] for index in pulse_indices]

    pi_times_us = np.full(len(normalized_steps), np.nan, dtype=float)
    pulse_pi_time_overrides = [normalized_steps[index]["pi_time_us"] for index in pulse_indices]
    missing_pulse_pi_override_indices = [local_index for local_index, value in enumerate(pulse_pi_time_overrides) if value is None]
    if missing_pulse_pi_override_indices:
        pulse_pi_times_from_inputs = None
        if per_pulse_pi_time_us is not None:
            provided_pi_times = np.asarray(per_pulse_pi_time_us, dtype=float)
            if provided_pi_times.shape == (len(normalized_steps),):
                pulse_pi_times_from_inputs = provided_pi_times[pulse_indices]
            elif provided_pi_times.shape == (len(pulse_indices),):
                pulse_pi_times_from_inputs = provided_pi_times
            else:
                raise ValueError("per_pulse_pi_time_us must have length equal to the step count or the driven-pulse count")
        elif transition_strengths_file is not None:
            if reference_pi_times_us is None or reference_strength_indices is None:
                raise ValueError("reference_pi_times_us and reference_strength_indices are required when using transition_strengths_file")
            strength_matrix = load_csv_matrix(str(transition_strengths_file))
            pi_time_matrix = compute_pi_time_matrix(strength_matrix, reference_pi_times_us, reference_strength_indices)
            pulse_pi_times_from_inputs = values_from_matrix(pulse_transitions, pi_time_matrix)
        else:
            pulse_pi_times_from_inputs = resolve_pi_times_us(pulse_couplings, default_pi_time_us, pi_time_by_coupling_us, None)

        for local_index, step_index in enumerate(pulse_indices):
            if normalized_steps[step_index]["pi_time_us"] is None:
                pi_times_us[step_index] = float(pulse_pi_times_from_inputs[local_index])
            else:
                pi_times_us[step_index] = float(normalized_steps[step_index]["pi_time_us"])
    else:
        for step_index in pulse_indices:
            pi_times_us[step_index] = float(normalized_steps[step_index]["pi_time_us"])

    sensitivities = np.full(len(normalized_steps), np.nan, dtype=float)
    if pulse_sensitivities_mhz_per_g is not None:
        provided_sensitivities = np.asarray(pulse_sensitivities_mhz_per_g, dtype=float)
        if provided_sensitivities.shape == (len(normalized_steps),):
            sensitivities[:] = provided_sensitivities
        elif provided_sensitivities.shape == (len(pulse_indices),):
            for local_index, step_index in enumerate(pulse_indices):
                sensitivities[step_index] = float(provided_sensitivities[local_index])
        else:
            raise ValueError("pulse_sensitivities_mhz_per_g must have length equal to the step count or the driven-pulse count")

    for step_index, step in enumerate(normalized_steps):
        if step["sensitivity_mhz_per_g"] is not None:
            sensitivities[step_index] = float(step["sensitivity_mhz_per_g"])

    unresolved_sensitivity_indices = [index for index, value in enumerate(sensitivities) if not np.isfinite(value)]
    if unresolved_sensitivity_indices:
        if sensitivity_matrix_file is None:
            raise ValueError("Wait steps require either sensitivity_matrix_file or a per-step sensitivity_mhz_per_g override")
        sensitivity_matrix = load_csv_matrix(str(sensitivity_matrix_file))
        resolved_sensitivities = values_from_matrix([transitions[index] for index in unresolved_sensitivity_indices], sensitivity_matrix)
        for step_index, value in zip(unresolved_sensitivity_indices, resolved_sensitivities):
            sensitivities[step_index] = float(value)

    durations_us = np.zeros(len(normalized_steps), dtype=float)
    if pulse_indices:
        durations_us[pulse_indices] = pulse_durations_us_from_thetas(thetas[pulse_indices], pi_times_us[pulse_indices])
    for step_index in wait_indices:
        durations_us[step_index] = float(normalized_steps[step_index]["duration_us"])
    timing_us = pulse_timing_us(durations_us)
    return {
        "dimension": dimension,
        "steps": normalized_steps,
        "step_kinds": step_kinds,
        "pulse_indices": pulse_indices,
        "wait_indices": wait_indices,
        "couplings": couplings,
        "thetas": thetas,
        "programmed_phases": phases,
        "transitions": transitions,
        "pi_times_us": pi_times_us,
        "pulse_sensitivities_mhz_per_g": sensitivities,
        "pulse_durations_us": durations_us,
        "durations_us": durations_us,
        "timing_us": timing_us,
    }


def run_sequence_simulation(prepared_sequence, initial_state, final_frame=None, apply_final_frame=False, simulation_mode="ideal", use_line_signal=True, line_model=None, line_reference="start", num_shots=200, rng_seed=None, calibration_scale=1.0, magnetic_scale=1.0, laser_scale=1.0, calibration_noise_hwhm_hz=15.0, calibration_noise_sigma_hz=0.0, magnetic_field_std_g=11.13e-6, laser_gaussian_T_s=2 * np.pi * 2206.8675860990334e-6, laser_lorentzian_L_s=2 * np.pi * 3353.155420102816e-6, voigt_x_min_hz=-1000.0, voigt_x_max_hz=1000.0, voigt_num_points=100000, rabi_frequency_std_fraction=0.037, line_signal_noise_std_g=0.0, target_state=None):
    simulation_mode = str(simulation_mode).strip().lower()
    initial_state = normalize_state(initial_state)
    deterministic = deterministic_sequence_result(
        prepared_sequence["couplings"],
        prepared_sequence["thetas"],
        prepared_sequence["programmed_phases"],
        initial_state=initial_state,
        pulse_sensitivities_mhz_per_g=prepared_sequence["pulse_sensitivities_mhz_per_g"],
        pi_times_us=prepared_sequence["pi_times_us"],
        durations_us=prepared_sequence["durations_us"],
        final_frame=final_frame,
        apply_final_frame=apply_final_frame,
        use_line_signal=use_line_signal,
        line_model=line_model,
        line_reference=line_reference,
        target_state=target_state,
    )
    if simulation_mode == "ideal":
        return {
            "mode": "ideal",
            "use_line_signal": bool(use_line_signal),
            "pi_times_us": deterministic["pi_times_us"],
            "pulse_durations_us": deterministic["pulse_durations_us"],
            "timing_us": deterministic["timing_us"],
            "line_phase_offsets_rad": deterministic["line_phase_offsets_rad"],
            "line_detunings_mhz": deterministic["line_detunings_mhz"],
            "line_reference_times_s": deterministic["line_reference_times_s"],
            "noiseless_state": deterministic["final_state"],
            "mean_populations": deterministic["populations"],
            "std_populations": np.zeros_like(deterministic["populations"]),
            "shot_populations": deterministic["populations"][None, :],
            "mean_fidelity": deterministic["fidelity"],
            "std_fidelity": 0.0 if np.isfinite(deterministic["fidelity"]) else np.nan,
            "deterministic_result": deterministic,
        }
    if simulation_mode != "monte_carlo":
        raise ValueError("simulation_mode must be 'ideal' or 'monte_carlo'")
    monte_carlo = monte_carlo_sequence_result(
        couplings=prepared_sequence["couplings"],
        thetas=prepared_sequence["thetas"],
        programmed_phases=prepared_sequence["programmed_phases"],
        initial_state=initial_state,
        pulse_sensitivities_mhz_per_g=prepared_sequence["pulse_sensitivities_mhz_per_g"],
        pi_times_us=prepared_sequence["pi_times_us"],
        durations_us=prepared_sequence["durations_us"],
        final_frame=final_frame,
        apply_final_frame=apply_final_frame,
        use_line_signal=use_line_signal,
        line_model=line_model,
        line_reference=line_reference,
        num_shots=num_shots,
        rng_seed=rng_seed,
        calibration_scale=calibration_scale,
        magnetic_scale=magnetic_scale,
        laser_scale=laser_scale,
        calibration_noise_hwhm_hz=calibration_noise_hwhm_hz,
        calibration_noise_sigma_hz=calibration_noise_sigma_hz,
        magnetic_field_std_g=magnetic_field_std_g,
        laser_gaussian_T_s=laser_gaussian_T_s,
        laser_lorentzian_L_s=laser_lorentzian_L_s,
        voigt_x_min_hz=voigt_x_min_hz,
        voigt_x_max_hz=voigt_x_max_hz,
        voigt_num_points=voigt_num_points,
        rabi_frequency_std_fraction=rabi_frequency_std_fraction,
        line_signal_noise_std_g=line_signal_noise_std_g,
        target_state=target_state,
    )
    monte_carlo["mode"] = "monte_carlo"
    monte_carlo["use_line_signal"] = bool(use_line_signal)
    monte_carlo["deterministic_result"] = deterministic
    return monte_carlo


def pulse_summary_records(prepared_sequence, result, mappings):
    line_phase_offsets = np.asarray(result["line_phase_offsets_rad"], dtype=float)
    line_detunings = np.asarray(result["line_detunings_mhz"], dtype=float)
    records = []
    for step_index, (step_kind, coupling, transition, theta, phase, pi_time_us, duration_us, timing, sensitivity, line_phase, line_detuning) in enumerate(
        zip(
            prepared_sequence["step_kinds"],
            prepared_sequence["couplings"],
            prepared_sequence["transitions"],
            prepared_sequence["thetas"],
            prepared_sequence["programmed_phases"],
            prepared_sequence["pi_times_us"],
            prepared_sequence["pulse_durations_us"],
            prepared_sequence["timing_us"],
            prepared_sequence["pulse_sensitivities_mhz_per_g"],
            line_phase_offsets,
            line_detunings,
        ),
        start=1,
    ):
        records.append(
            {
                "step": step_index,
                "kind": step_kind,
                "coupling": format_coupling_label(coupling, mappings),
                "transition": format_transition_label(transition),
                "theta_rad": float(theta),
                "theta_over_pi": float(theta / np.pi),
                "phase_rad": float(phase),
                "pi_time_us": None if not np.isfinite(pi_time_us) else float(pi_time_us),
                "duration_us": float(duration_us),
                "t_start_us": float(timing[0]),
                "t_end_us": float(timing[1]),
                "sensitivity_MHz_per_G": float(sensitivity),
                "line_phase_rad": float(line_phase),
                "line_detuning_MHz": float(line_detuning),
            }
        )
    return records


def population_summary_records(result, mappings):
    mean_populations = np.asarray(result["mean_populations"], dtype=float)
    std_populations = np.asarray(result["std_populations"], dtype=float)
    records = []
    for state_index, (mean_value, std_value) in enumerate(zip(mean_populations, std_populations)):
        records.append(
            {
                "state_index": state_index,
                "state_label": format_state_label(state_index, mappings),
                "mean_population": float(mean_value),
                "std_population": float(std_value),
            }
        )
    return records


def summarize_simulation_result(result, prepared_sequence, mappings):
    print(f"mode = {result['mode']}")
    print(f"use_line_signal = {result['use_line_signal']}")
    print(f"step_count = {len(prepared_sequence['couplings'])}")
    print(f"driven_pulse_count = {len(prepared_sequence.get('pulse_indices', []))}")
    print(f"wait_count = {len(prepared_sequence.get('wait_indices', []))}")
    print(f"total_sequence_time_us = {np.sum(prepared_sequence['durations_us']):.6f}")
    print(f"line_phase_offsets_rad = {result['line_phase_offsets_rad']}")
    print(f"line_detunings_mhz = {result['line_detunings_mhz']}")
    if np.isfinite(result["mean_fidelity"]):
        print(f"mean_fidelity = {result['mean_fidelity']}")
        if result["mode"] == "monte_carlo":
            print(f"std_fidelity = {result['std_fidelity']}")
    display_records(population_summary_records(result, mappings))


def plot_populations(result, mappings, title_prefix="Sequence simulation"):
    mean_populations = 100.0 * np.asarray(result["mean_populations"], dtype=float)
    std_populations = 100.0 * np.asarray(result["std_populations"], dtype=float)
    state_labels = [format_state_label(index, mappings) for index in range(len(mean_populations))]
    x_values = np.arange(len(mean_populations))
    fig_width = max(8.0, 0.8 * len(mean_populations))
    fig, ax = plt.subplots(figsize=(fig_width, 4.8))
    ax.bar(x_values, mean_populations, color="#1f4e79", alpha=0.85, label="Mean population")
    if np.any(std_populations > 0):
        ax.errorbar(x_values, mean_populations, yerr=std_populations, fmt="none", ecolor="black", elinewidth=1.2, capsize=4)
    if result["mode"] == "monte_carlo":
        deterministic = 100.0 * np.asarray(result["deterministic_result"]["populations"], dtype=float)
        ax.plot(x_values, deterministic, "kx", markersize=8, label="Deterministic baseline")
    ax.set_xticks(x_values)
    ax.set_xticklabels(state_labels, rotation=45, ha="right")
    ax.set_ylabel("Population (%)")
    ax.set_xlabel("State")
    ax.set_ylim(0.0, 100.0)
    ax.set_title(f"{title_prefix} ({result['mode']})")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend()
    plt.tight_layout()
    plt.show()
