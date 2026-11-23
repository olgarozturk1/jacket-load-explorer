"""Save a ten-second AVI from the generated load-cycle data.

Run build_viewer.py first, then python export_avi.py.
AVI has a fixed camera; the HTML viewer remains rotatable.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, writers


def export_avi():
    root = Path(__file__).resolve().parent
    data_file = root / "generated" / "load_cycle.npz"
    if not data_file.exists():
        raise FileNotFoundError("Run build_viewer.py first to calculate the animation data.")

    with np.load(data_file) as data:
        headings = data["headings"]
        phases = data["phases"]
        loads = data["load_MN"]
        period = float(data["period"])

    n = len(headings)
    angles = np.append(headings, 360.0)
    closed = np.zeros((len(phases), n + 1, n + 1))
    closed[:, :-1, :-1] = loads
    closed[:, -1, :-1] = loads[:, 0, :]
    closed[:, :-1, -1] = loads[:, :, 0]
    closed[:, -1, -1] = loads[:, 0, 0]
    waves, currents = np.meshgrid(angles, angles, indexing="ij")
    ceiling = max(1.0, np.ceil(loads.max()))
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)

    def draw(frame):
        ax.clear()
        surface = ax.plot_surface(
            currents, waves, closed[frame], cmap="viridis",
            vmin=0, vmax=ceiling, rstride=1, cstride=1,
            linewidth=0.1, edgecolor=(0, 0, 0, 0.12), zorder=1,
        )
        i, j = np.unravel_index(np.argmax(loads[frame]), (n, n))
        ax.scatter(
            headings[j], headings[i], loads[frame, i, j],
            color="red", edgecolors="white", s=65,
            depthshade=False, zorder=10,
        )
        ax.set_title(
            f"Instantaneous horizontal load | Phase {phases[frame]:.1f} deg | "
            f"t = {phases[frame] / 360 * period:.2f} s", pad=22,
        )
        ax.set_xlabel("Current direction (degrees)", labelpad=10)
        ax.set_ylabel("Wave direction (degrees)", labelpad=10)
        ax.set_zlabel("Horizontal resultant (MN)", labelpad=10)
        ax.set_xlim(0, 360)
        ax.set_ylim(0, 360)
        ax.set_zlim(0, ceiling)
        ax.set_xticks([0, 90, 180, 270, 360])
        ax.set_yticks([0, 90, 180, 270, 360])
        ax.view_init(elev=28, azim=-125)
        return surface

    surface = draw(0)
    fig.colorbar(surface, ax=ax, shrink=0.6, pad=0.12, label="Load (MN)")
    fig.text(
        0.5, 0.04, "Red dot: largest sampled load at this phase. Phase is measured at the origin.",
        ha="center", fontsize=9,
    )
    movie = FuncAnimation(fig, draw, frames=len(phases), blit=False)
    if not writers.is_available("ffmpeg"):
        import imageio_ffmpeg
        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    writer = FFMpegWriter(
        fps=len(phases) / period, codec="mjpeg",
        extra_args=["-q:v", "3", "-pix_fmt", "yuvj420p"],
    )
    output = root / "generated" / "jacket_load_cycle.avi"
    print("Saving AVI — please wait...")
    movie.save(str(output), writer=writer, dpi=100)
    plt.close(fig)
    print(f"Saved: {output}")


if __name__ == "__main__":
    export_avi()
