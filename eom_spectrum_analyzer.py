"""
EOM (Electro-Optic Modulator) Spectrum Analyzer
Simulates the optical spectrum when multiple RF tones are applied to an EOM.
Models phase and amplitude modulation with sideband generation.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv  # Bessel function of the first kind
from scipy.fft import fft, fftfreq

class EOMSimulator:
    """
    Simulates an Electro-Optic Modulator with multiple RF driving signals.
    """
    
    def __init__(self, carrier_wavelength=780e-9, sampling_rate=100e9, 
                 duration=10e-6, num_points=100000):
        """
        Initialize the EOM simulator.
        
        Parameters:
        -----------
        carrier_wavelength : float
            Optical carrier wavelength in meters (default: 780 nm for Rb)
        sampling_rate : float
            Sampling rate in Hz (default: 100 GHz)
        duration : float
            Signal duration in seconds (default: 10 microseconds)
        num_points : int
            Number of time points
        """
        self.carrier_wavelength = carrier_wavelength
        self.carrier_frequency = 3e8 / carrier_wavelength  # c / λ
        self.sampling_rate = sampling_rate
        self.duration = duration
        self.num_points = num_points
        self.time = np.linspace(0, duration, num_points)
        
        # RF modulation signals
        self.rf_signals = []
        
    def add_rf_tone(self, frequency, modulation_depth, phase=0, name=None):
        """
        Add an RF tone to modulate the EOM.
        
        Parameters:
        -----------
        frequency : float
            RF frequency in Hz
        modulation_depth : float
            Modulation depth β (dimensionless, typically 0-5)
            β = π * V_rf / V_π for phase modulators
        phase : float
            Phase in degrees (default: 0)
        name : str
            Optional name for the RF tone
        """
        rf_dict = {
            'frequency': frequency,
            'modulation_depth': modulation_depth,
            'phase': phase,
            'name': name or f'RF Tone {len(self.rf_signals) + 1}'
        }
        self.rf_signals.append(rf_dict)
        
    def clear_rf_tones(self):
        """Clear all RF tones."""
        self.rf_signals = []
        
    def generate_modulation_signal(self):
        """
        Generate the total modulation signal from all RF tones.
        
        Returns:
        --------
        numpy.ndarray
            Time-varying modulation signal β(t)
        """
        modulation = np.zeros_like(self.time)
        
        for rf in self.rf_signals:
            freq = rf['frequency']
            beta = rf['modulation_depth']
            phase_rad = np.deg2rad(rf['phase'])
            
            # Add this RF tone to the modulation
            modulation += beta * np.cos(2 * np.pi * freq * self.time + phase_rad)
        
        return modulation
    
    def generate_optical_field(self, modulation_type='phase'):
        """
        Generate the modulated optical field.
        
        Parameters:
        -----------
        modulation_type : str
            'phase' for phase modulation or 'amplitude' for amplitude modulation
            
        Returns:
        --------
        numpy.ndarray
            Complex electric field E(t)
        """
        # Get total modulation signal
        modulation = self.generate_modulation_signal()
        
        # Carrier
        carrier_phase = 2 * np.pi * self.carrier_frequency * self.time
        
        if modulation_type == 'phase':
            # Phase modulation: E(t) = E₀ exp[i(ω_c*t + β(t))]
            field = np.exp(1j * (carrier_phase + modulation))
        else:  # amplitude modulation
            # Amplitude modulation: E(t) = E₀[1 + m(t)] exp[i*ω_c*t]
            field = (1 + modulation) * np.exp(1j * carrier_phase)
        
        return field
    
    def compute_optical_spectrum(self, modulation_type='phase', max_sidebands=10):
        """
        Compute the optical spectrum using analytical Bessel expansion.
        More accurate than FFT for narrow-line analysis.
        
        Parameters:
        -----------
        modulation_type : str
            'phase' or 'amplitude' modulation
        max_sidebands : int
            Maximum number of sidebands to compute for each tone
            
        Returns:
        --------
        frequencies : numpy.ndarray
            Frequency array relative to carrier (Hz)
        spectrum : numpy.ndarray
            Spectrum amplitudes (normalized)
        """
        # For single tone, use Bessel function expansion
        if len(self.rf_signals) == 1:
            rf = self.rf_signals[0]
            freq_rf = rf['frequency']
            beta = rf['modulation_depth']
            
            # Calculate sidebands using Bessel functions
            orders = np.arange(-max_sidebands, max_sidebands + 1)
            frequencies = orders * freq_rf
            
            if modulation_type == 'phase':
                # Phase modulation: J_n(β) for nth sideband
                spectrum = np.abs(jv(orders, beta))
            else:
                # Amplitude modulation: different pattern
                spectrum = np.zeros_like(frequencies, dtype=float)
                spectrum[orders == 0] = 1.0
                spectrum[orders == 1] = beta / 2
                spectrum[orders == -1] = beta / 2
                
        else:
            # Multiple tones: use FFT method
            field = self.generate_optical_field(modulation_type)
            
            # Compute intensity (|E|²)
            intensity = np.abs(field)**2
            
            # FFT of intensity to get beat notes
            spectrum_fft = fft(intensity - np.mean(intensity))
            freqs_fft = fftfreq(len(intensity), 1/self.sampling_rate)
            
            # Only positive frequencies
            positive_idx = freqs_fft >= 0
            frequencies = freqs_fft[positive_idx]
            spectrum = np.abs(spectrum_fft[positive_idx])
            
            # Normalize
            spectrum = spectrum / np.max(spectrum) if np.max(spectrum) > 0 else spectrum
            
        return frequencies, spectrum
    
    def compute_fft_spectrum(self, modulation_type='phase'):
        """
        Compute optical spectrum using FFT (good for time-domain analysis).
        
        Parameters:
        -----------
        modulation_type : str
            'phase' or 'amplitude' modulation
            
        Returns:
        --------
        frequencies : numpy.ndarray
            Frequency array (Hz)
        spectrum : numpy.ndarray
            Power spectrum (normalized, dB scale)
        """
        # Generate field
        field = self.generate_optical_field(modulation_type)
        
        # FFT
        spectrum_fft = fft(field)
        freqs = fftfreq(len(field), 1/self.sampling_rate)
        
        # Shift to center carrier
        center_idx = len(freqs) // 2
        freqs_shifted = freqs - self.carrier_frequency
        
        # Power spectrum
        power = np.abs(spectrum_fft)**2
        
        # Only show around carrier (±50 GHz window)
        window = 50e9
        mask = np.abs(freqs_shifted) < window
        
        frequencies = freqs_shifted[mask]
        spectrum = power[mask]
        
        # Normalize and convert to dB
        spectrum = spectrum / np.max(spectrum) if np.max(spectrum) > 0 else spectrum
        spectrum_db = 10 * np.log10(spectrum + 1e-12)
        
        return frequencies, spectrum_db
    
    def plot_spectrum(self, modulation_type='phase', method='analytical', max_sidebands=10, 
                     figsize=(14, 10)):
        """
        Plot the optical spectrum and modulation signals.
        
        Parameters:
        -----------
        modulation_type : str
            'phase' or 'amplitude' modulation
        method : str
            'analytical' (Bessel) or 'fft' for spectrum calculation
        max_sidebands : int
            Maximum sidebands for analytical method
        figsize : tuple
            Figure size
        """
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(3, 1, hspace=0.3)
        
        # Plot 1: RF modulation signals
        ax1 = fig.add_subplot(gs[0])
        modulation = self.generate_modulation_signal()
        time_ns = self.time * 1e9
        
        ax1.plot(time_ns, modulation, color='blue', linewidth=1.5)
        ax1.set_xlabel('Time (ns)')
        ax1.set_ylabel('Modulation β(t)')
        ax1.set_title('RF Modulation Signal (Combined)')
        ax1.grid(True, alpha=0.3)
        
        # Add text with RF tone info
        info_text = f"RF Tones: {len(self.rf_signals)}\n"
        for rf in self.rf_signals:
            info_text += f"{rf['name']}: {rf['frequency']/1e6:.1f} MHz, β={rf['modulation_depth']:.2f}\n"
        ax1.text(0.02, 0.98, info_text, transform=ax1.transAxes, 
                verticalalignment='top', fontsize=8, bbox=dict(boxstyle='round', 
                facecolor='wheat', alpha=0.5))
        
        # Plot 2: Optical field (real part, zoomed in)
        ax2 = fig.add_subplot(gs[1])
        field = self.generate_optical_field(modulation_type)
        
        # Show only first 100 ns for visibility
        plot_duration = min(100e-9, self.duration)
        plot_points = int(plot_duration * self.sampling_rate)
        time_plot = self.time[:plot_points] * 1e9
        
        ax2.plot(time_plot, np.real(field[:plot_points]), color='red', linewidth=0.5)
        ax2.set_xlabel('Time (ns)')
        ax2.set_ylabel('E(t) [Re]')
        ax2.set_title(f'Modulated Optical Field ({modulation_type.capitalize()} Modulation)')
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Optical spectrum
        ax3 = fig.add_subplot(gs[2])
        
        if method == 'analytical':
            frequencies, spectrum = self.compute_optical_spectrum(modulation_type, max_sidebands)
            
            # Stem plot for discrete sidebands
            markerline, stemlines, baseline = ax3.stem(frequencies/1e6, spectrum, 
                                                        linefmt='red', markerfmt='ro',
                                                        basefmt='k-')
            markerline.set_markersize(6)
            stemlines.set_linewidth(1.5)
            
            ax3.set_ylabel('Amplitude (normalized)')
            ax3.set_ylim([0, max(1.1, np.max(spectrum) * 1.1)])
            
        else:  # FFT method
            frequencies, spectrum_db = self.compute_fft_spectrum(modulation_type)
            ax3.plot(frequencies/1e6, spectrum_db, color='red', linewidth=1.5)
            ax3.set_ylabel('Power (dB)')
            ax3.set_ylim([-80, 5])
        
        ax3.set_xlabel('Frequency Relative to Carrier (MHz)')
        ax3.set_title(f'Optical Spectrum (Method: {method})')
        ax3.grid(True, alpha=0.3)
        ax3.axvline(0, color='blue', linestyle='--', alpha=0.5, linewidth=2, label='Carrier')
        
        # Mark RF tone positions
        for rf in self.rf_signals:
            freq_mhz = rf['frequency'] / 1e6
            ax3.axvline(freq_mhz, color='green', linestyle=':', alpha=0.5)
            ax3.axvline(-freq_mhz, color='green', linestyle=':', alpha=0.5)
        
        ax3.legend()
        
        plt.tight_layout()
        return fig


def run_from_lists(rf_frequencies, modulation_depths, rf_phases=None, rf_names=None,
                   carrier_wavelength=780e-9, modulation_type='phase', 
                   sampling_rate=100e9, duration=10e-6, num_points=100000,
                   plot_method='analytical', max_sidebands=10):
    """
    Run EOM simulator with lists of RF parameters.
    
    Parameters:
    -----------
    rf_frequencies : list of float
        List of RF frequencies in Hz
    modulation_depths : list of float
        List of modulation depths β (dimensionless)
    rf_phases : list of float, optional
        List of phases in degrees (default: all zeros)
    rf_names : list of str, optional
        List of RF tone names
    carrier_wavelength : float
        Optical carrier wavelength in meters (default: 780 nm)
    modulation_type : str
        'phase' or 'amplitude' modulation
    sampling_rate : float
        Sampling rate in Hz
    duration : float
        Signal duration in seconds
    num_points : int
        Number of time points
    plot_method : str
        'analytical' or 'fft'
    max_sidebands : int
        Maximum sidebands for analytical method
        
    Returns:
    --------
    fig : matplotlib.figure.Figure
        Generated figure
    eom : EOMSimulator
        EOM simulator object
    """
    # Validate inputs
    if not (len(rf_frequencies) == len(modulation_depths)):
        raise ValueError("rf_frequencies and modulation_depths must have same length")
    
    if rf_phases is None:
        rf_phases = [0] * len(rf_frequencies)
    
    # Create EOM
    eom = EOMSimulator(
        carrier_wavelength=carrier_wavelength,
        sampling_rate=sampling_rate,
        duration=duration,
        num_points=num_points
    )
    
    # Add RF tones
    for i, (freq, beta, phase) in enumerate(zip(rf_frequencies, modulation_depths, rf_phases)):
        name = rf_names[i] if rf_names and i < len(rf_names) else f"RF{i+1}"
        eom.add_rf_tone(frequency=freq, modulation_depth=beta, phase=phase, name=name)
    
    # Print info
    print(f"\nEOM Configuration:")
    print(f"Carrier wavelength: {carrier_wavelength*1e9:.1f} nm")
    print(f"Carrier frequency: {eom.carrier_frequency/1e12:.3f} THz")
    print(f"Modulation type: {modulation_type}")
    print(f"\nRF Tones ({len(rf_frequencies)}):")
    for i, (freq, beta) in enumerate(zip(rf_frequencies, modulation_depths)):
        print(f"  {i+1}. Frequency: {freq/1e6:.3f} MHz, β: {beta:.3f}")
    
    # Plot
    fig = eom.plot_spectrum(modulation_type=modulation_type, method=plot_method, 
                           max_sidebands=max_sidebands)
    
    return fig, eom


if __name__ == "__main__":
    # ============================================================================
    # CONFIGURATION: Edit these parameters for your EOM
    # ============================================================================
    
    # Optical carrier
    carrier_wavelength = 780e-9  # 780 nm (Rubidium D2 line)
    modulation_type = 'phase'     # 'phase' or 'amplitude'
    
    # RF tones applied to EOM
    rf_frequencies = [10e6, 20e6]        # RF frequencies in Hz (e.g., 10 MHz, 20 MHz)
    modulation_depths = [1.5, 0.8]       # Modulation depths β (β = π*V_rf/V_π)
    rf_phases = [0, 0]                   # Phases in degrees
    rf_names = ['Tone 1', 'Tone 2']      # Optional names
    
    # Simulation parameters
    sampling_rate = 100e9    # 100 GHz sampling
    duration = 10e-6         # 10 microseconds
    num_points = 100000      # Time points
    
    # Plotting
    plot_method = 'analytical'  # 'analytical' (Bessel) or 'fft'
    max_sidebands = 10          # Max sidebands for analytical method
    
    # ============================================================================
    # RUN
    # ============================================================================
    
    print("="*70)
    print("EOM Spectrum Analyzer")
    print("="*70)
    
    fig, eom = run_from_lists(
        rf_frequencies=rf_frequencies,
        modulation_depths=modulation_depths,
        rf_phases=rf_phases,
        rf_names=rf_names,
        carrier_wavelength=carrier_wavelength,
        modulation_type=modulation_type,
        sampling_rate=sampling_rate,
        duration=duration,
        num_points=num_points,
        plot_method=plot_method,
        max_sidebands=max_sidebands
    )
    
    plt.show()
