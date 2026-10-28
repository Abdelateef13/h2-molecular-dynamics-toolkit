"""Research-friendly plots and output files."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .analysis import DiffusionFit
from .trajectory import Trajectory


def plot_analysis(trajectory: Trajectory, lag_times: np.ndarray, msd: np.ndarray,
                  fit: DiffusionFit, radii: np.ndarray, rdf: np.ndarray,
                  output_dir: str | Path) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 125})
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    ax = axs[0]
    ax.plot(lag_times, msd, color="#067A96", lw=2.3, label="Time-origin averaged MSD")
    mask = (lag_times >= fit.start_time) & (lag_times <= fit.end_time)
    ax.plot(lag_times[mask], fit.slope * lag_times[mask] + fit.intercept,
            "--", color="#E97B28", lw=2.3, label=f"Linear fit ($R^2$={fit.r_squared:.3f})")
    ax.set(xlabel=f"Lag time ({trajectory.time_unit})",
           ylabel=f"MSD ({trajectory.length_unit}$^2$)", title="H₂ translational motion")
    ax.legend(fontsize=9, frameon=False)
    ax.grid(alpha=0.17)
    ax = axs[1]
    ax.plot(radii, rdf, color="#264D87", lw=2.0)
    ax.axhline(1, ls="--", color="#888888", alpha=0.7, label="Ideal-gas baseline")
    ax.set(xlabel=f"Distance ({trajectory.length_unit})", ylabel="g(r)",
           title="H₂–H₂ radial distribution")
    ax.grid(alpha=0.17)
    ax.legend(fontsize=9, frameon=False)
    fig.suptitle("H₂-only molecular dynamics | structural & transport analysis", fontsize=14)
    fig.savefig(out / "h2_analysis.png", dpi=190, bbox_inches="tight")
    plt.close(fig)

    if trajectory.potential_energy is not None and trajectory.kinetic_energy is not None:
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
        t = trajectory.times
        total = trajectory.potential_energy + trajectory.kinetic_energy
        axs[0].plot(t, total, color="#067A96", lw=1.7, label="Total (NVE)")
        axs[0].plot(t, trajectory.kinetic_energy, color="#E97B28", lw=1.4,
                    label="Kinetic")
        axs[0].plot(t, trajectory.potential_energy, color="#264D87", lw=1.4,
                    label="Potential")
        axs[0].set(xlabel=f"Production time ({trajectory.time_unit})",
                   ylabel="Energy ($\\epsilon$)", title="Energy conservation")
        axs[0].legend(frameon=False, fontsize=9)
        axs[0].grid(alpha=.17)
        if trajectory.temperature is not None:
            axs[1].plot(t, trajectory.temperature, color="#264D87", lw=1.5)
            axs[1].axhline(trajectory.temperature.mean(), ls="--", color="#E97B28",
                           label="Production mean")
            axs[1].set(xlabel=f"Production time ({trajectory.time_unit})",
                       ylabel="Temperature ($\\epsilon/k_B$)", title="NVE temperature")
            axs[1].grid(alpha=.17)
            axs[1].legend(frameon=False, fontsize=9)
        fig.suptitle("Simulation diagnostics | illustrative reduced-unit LJ model", fontsize=13)
        fig.savefig(out / "h2_diagnostics.png", dpi=190, bbox_inches="tight")
        plt.close(fig)
