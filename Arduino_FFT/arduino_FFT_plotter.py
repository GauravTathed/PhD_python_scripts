import time
import serial
import numpy as np
import matplotlib.pyplot as plt

PORT = 'COM3'
BAUD = 2000000
CHUNK = 2048
BUF_LEN = 16384

ADC_MAX = 4095.0
VREF = 3.3
R_LOAD = 50.0

ser = serial.Serial(PORT, BAUD, timeout=1)

plt.ion()
fig, (ax_time, ax_fft) = plt.subplots(2, 1)

x_time = np.arange(BUF_LEN)
y_time = np.zeros(BUF_LEN, dtype=np.float32)
line_time, = ax_time.plot(x_time, y_time)
ax_time.set_xlim(0, BUF_LEN)
ax_time.set_ylim(0, VREF)
ax_time.set_xlabel('Sample index')
ax_time.set_ylabel('Voltage (V)')

text_vpp = ax_time.text(0.98, 0.02, "Ripple Vpp(ac): 0.0000 V",
                        transform=ax_time.transAxes,
                        ha='right', va='bottom')

line_fft, = ax_fft.plot([], [])
ax_fft.set_xlim(0, 360)
ax_fft.set_xlabel('Frequency (Hz)')
ax_fft.set_ylabel('Magnitude (dBm)')

text_markers = ax_fft.text(0.98, 0.02, "",
                           transform=ax_fft.transAxes,
                           ha='right', va='bottom')

plt.show(block=False)

prev_t = time.perf_counter()
fs_est = None

while True:
    raw = ser.read(2 * CHUNK)
    if len(raw) != 2 * CHUNK:
        plt.pause(0.001)
        continue

    now = time.perf_counter()
    dt = now - prev_t
    prev_t = now
    if dt > 0:
        fs_est = CHUNK / dt

    adc = np.frombuffer(raw, dtype='<u2').astype(np.float32)
    volts = (adc / ADC_MAX) * VREF

    y_time = np.roll(y_time, -CHUNK)
    y_time[-CHUNK:] = volts
    line_time.set_ydata(y_time)

    y_ac = y_time - np.mean(y_time)
    vpp = float(y_ac.max() - y_ac.min())
    text_vpp.set_text(f"Ripple Vpp(ac): {vpp:.4f} V")

    if fs_est is not None and fs_est > 0:
        y = y_time - np.mean(y_time)
        fft_vals = np.fft.rfft(y)
        N = BUF_LEN
        mag = np.abs(fft_vals) / N
        if len(mag) > 1:
            mag[1:] *= 2.0
        vrms_bins = mag / np.sqrt(2.0)
        P_watt = (vrms_bins ** 2) / R_LOAD
        P_watt = np.maximum(P_watt, 1e-15)
        dBm = 10.0 * np.log10(P_watt / 1e-3)

        freqs = np.fft.rfftfreq(N, d=1.0 / fs_est)
        mask = freqs <= 360.0
        f_plot = freqs[mask]
        dBm_plot = dBm[mask]

        if len(f_plot) > 0:
            line_fft.set_data(f_plot, dBm_plot)
            ymin = float(np.min(dBm_plot)) - 3.0
            ymax = float(np.max(dBm_plot)) + 3.0
            ax_fft.set_xlim(0, 360)
            ax_fft.set_ylim(ymin, ymax)

            def pick_freq(target):
                idx = int(np.argmin(np.abs(f_plot - target)))
                if 0 <= idx < len(dBm_plot):
                    return f_plot[idx], dBm_plot[idx]
                return None, None

            f60, d60 = pick_freq(60.0)
            f180, d180 = pick_freq(180.0)

            s60 = "N/A"
            s180 = "N/A"
            if d60 is not None:
                s60 = f"{f60:.1f} Hz: {d60:.1f} dBm"
            if d180 is not None:
                s180 = f"{f180:.1f} Hz: {d180:.1f} dBm"

            text_markers.set_text(f"60 Hz: {s60}\n180 Hz: {s180}")

    fig.canvas.draw_idle()
    plt.pause(0.001)
