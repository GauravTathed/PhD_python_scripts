import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import imageio.v2 as imageio

def line_signal_model(t_us, B0, A60, phi60, A180, phi180):
    t_s = t_us * 1e-6
    return (
        B0
        + A60 * np.cos(2 * np.pi * 60 * t_s + phi60)
        + A180 * np.cos(2 * np.pi * 180 * t_s + phi180)
    )

def save_simple_gif(times_s, pops, B_func, out_gif="test.gif", fps=30, dpi=120):
    t_us = times_s / 1e-6
    B_mG = np.array([B_func(t) for t in times_s]) * 1e3

    fig, ax = plt.subplots(figsize=(9, 5), dpi=dpi)

    (l0,) = ax.plot([], [], lw=2, label="|0⟩")
    (l1,) = ax.plot([], [], lw=2, label="|1⟩")
    (l2,) = ax.plot([], [], lw=2, label="|2⟩")

    ax.set_xlim(t_us[0], t_us[-1])
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel("Population")
    ax.grid(True)

    axB = ax.twinx()
    (lb,) = axB.plot([], [], "--", lw=1.8, label="B(t)")
    axB.set_ylabel("Magnetic field (mG)")
    span = np.ptp(B_mG)
    pad = 0.05 * span if span > 0 else 0.1
    axB.set_ylim(np.min(B_mG) - pad, np.max(B_mG) + pad)

    h1, lab1 = ax.get_legend_handles_labels()
    h2, lab2 = axB.get_legend_handles_labels()
    ax.legend(h1 + h2, lab1 + lab2, loc="upper right")

    frames = []

    def update(k):
        l0.set_data(t_us[:k+1], pops[:k+1, 0])
        l1.set_data(t_us[:k+1], pops[:k+1, 1])
        l2.set_data(t_us[:k+1], pops[:k+1, 2])
        lb.set_data(t_us[:k+1], B_mG[:k+1])
        ax.set_title(f"t = {t_us[k]:.2f} µs")

    for k in range(len(times_s)):
        update(k)
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba())
        frames.append(rgba[..., :3].copy())

    imageio.mimsave(out_gif, frames, fps=fps)
    plt.close(fig)

if __name__ == "__main__":
    # fake “populations” just to test saving
    times_s = np.linspace(0, 1e-3, 400)
    p0 = 0.5 * (1 + np.cos(2 * np.pi * 8e3 * times_s))
    p1 = 1 - p0
    p2 = np.zeros_like(p0)
    pops = np.column_stack([p0, p1, p2])

    # your fitted line-synchronous B field (mG -> converted inside B_func to G)
    B0_mG = 0.3076
    A60_mG = -0.2973
    phi60 = 0.8148
    A180_mG = -0.0969
    phi180 = -0.4964

    B_func = lambda t_s: line_signal_model(
        t_s / 1e-6, B0_mG, A60_mG, phi60, A180_mG, phi180
    ) * 1e-3  # returns G

    save_simple_gif(times_s, pops, B_func, out_gif="simple_test.gif", fps=30, dpi=120)
    print("Saved: simple_test.gif")
