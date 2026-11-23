"""
PART 2(a): HORIZONTAL FORCE ON THE WHOLE JACKET

Numerical methods:
    - np.linspace: sample positions, phases and directions
    - np.interp: locate the wave number and water-surface intersection
    - Trapezoidal rule: integrate forces
    - Loops and np.argmax: search for the maximum
    - Grid refinement: check numerical accuracy


Physical assumptions:
    - Independent wave and current headings
    - Linear wave-current theory
    - First-order material acceleration
    - Both drag and inertia use the instantaneous submerged length
    - Linear velocity profiles continued above SWL
    - Member-centreline wetting
    - Columns: vector horizontal drag
    - Beams/diagonals: independent normal-component drag
    - No axial drag, shielding or diffraction

Heading = direction of travel, measured from +x towards +y.
Phase = omega*t at the origin.
"""

import numpy as np
import matplotlib.pyplot as plt


# =====================================================================
# 1. INPUTS
# =====================================================================

DEPTH = 40.0
AMPLITUDE = 6.0
PERIOD = 10.0
CURRENT_SPEED = 1.6

RHO = 1025.0
GRAVITY = 9.81
CD = 1.2
CM = 0.8

OMEGA = 2.0 * np.pi / PERIOD

# Sample wave numbers for interpolation.
K_SAMPLES = np.linspace(0.0, 0.12, 12001)


# =====================================================================
# 2. GEOMETRY
# =====================================================================

GROUPS = [
    "Columns",           # group 0
    "Horizontal beams",  # group 1
    "x-z diagonals",     # group 2
    "y-z diagonals",     # group 3
]

# Each row:
# group, diameter, start coordinates, end coordinates
#
# Coordinates are (x, y, z), in metres.
# Columns and diagonals start at their LOWER end.

MEMBERS = [
    # Columns
    (0, 4.0, (0, 0, -40), (0, 0, 12)),
    (0, 4.0, (0, 30, -40), (0, 30, 12)),
    (0, 4.0, (30, 0, -40), (30, 0, 12)),
    (0, 4.0, (30, 30, -40), (30, 30, 12)),

    # Lower beams
    (1, 2.0, (0, 0, -35), (30, 0, -35)),
    (1, 2.0, (0, 30, -35), (30, 30, -35)),
    (1, 2.0, (0, 0, -35), (0, 30, -35)),
    (1, 2.0, (30, 0, -35), (30, 30, -35)),

    # Upper beams
    (1, 2.0, (0, 0, 7), (30, 0, 7)),
    (1, 2.0, (0, 30, 7), (30, 30, 7)),
    (1, 2.0, (0, 0, 7), (0, 30, 7)),
    (1, 2.0, (30, 0, 7), (30, 30, 7)),

    # x-z diagonals
    (2, 1.0, (0, 0, -35), (30, 0, 7)),
    (2, 1.0, (0, 30, -35), (30, 30, 7)),

    # y-z diagonals
    (3, 1.0, (0, 30, -35), (0, 0, 7)),
    (3, 1.0, (30, 30, -35), (30, 0, 7)),
]


# =====================================================================
# 3. FIND k AND sigma USING SAMPLING AND INTERPOLATION
# =====================================================================

def wave_parameters(wave_heading, current_heading):
    wave_angle = np.deg2rad(wave_heading)
    current_angle = np.deg2rad(current_heading)

    wave_direction = np.array([
        np.cos(wave_angle),
        np.sin(wave_angle),
    ])

    current_velocity = CURRENT_SPEED * np.array([
        np.cos(current_angle),
        np.sin(current_angle),
    ])

    # Current component along the wave direction.
    current_along_wave = np.sum(
        current_velocity * wave_direction
    )

    # Frequency corresponding to each sampled k:
    #
    # omega = sqrt(g*k*tanh(k*d)) + k*current_along_wave
    frequency = np.sqrt(
        GRAVITY * K_SAMPLES * np.tanh(K_SAMPLES * DEPTH)
    )

    frequency = (
        frequency + current_along_wave * K_SAMPLES
    )

    # np.interp requires increasing sample coordinates.
    if (
        np.any(np.diff(frequency) <= 0)
        or not frequency[0] <= OMEGA <= frequency[-1]
    ):
        raise ValueError(
            "The k sampling interval must be reconsidered for these inputs."
        )

    # Find the k corresponding to our known OMEGA.
    k = np.interp(
        OMEGA,
        frequency,
        K_SAMPLES,
    )

    sigma = OMEGA - k * current_along_wave

    return wave_direction, current_velocity, k, sigma


# =====================================================================
# 4. FIND THE SUBMERGED LENGTH UP TO eta
# =====================================================================

def submerged_fraction(
    start_z,
    rise,
    start_phase,
    phase_change,
    fractions,
):
    """
    Find the fraction of the member below the local water surface.

    Along the member:
        z = start_z + rise*fraction
        eta = amplitude*sin(start_phase - phase_change*fraction)

    The intersection occurs where z - eta = 0.
    """

    if rise == 0.0:
        # For this geometry:
        # lower beams are always wet;
        # upper beam centrelines are always dry.
        if start_z < -AMPLITUDE:
            return np.ones_like(start_phase)

        if start_z >= AMPLITUDE:
            return np.zeros_like(start_phase)

        raise ValueError(
            "A beam inside the wave envelope needs different wetting logic."
        )

    z = start_z + rise * fractions

    wet_fraction = np.zeros_like(start_phase)

    for i in range(len(start_phase)):
        eta = AMPLITUDE * np.sin(
            start_phase[i] - phase_change * fractions
        )

        gap = z - eta

        # For these rising columns/diagonals, gap increases along
        # the member. Interpolate to the point where gap = 0.
        wet_fraction[i] = np.interp(
            0.0,
            gap,
            fractions,
        )

    return wet_fraction


# =====================================================================
# 5. TRAPEZOIDAL INTEGRATION
# =====================================================================

def integrate_force(force_per_metre, distance):
    """
    Chapter 4 trapezoidal rule:

        integral ≈ sum[
            0.5*(f[i] + f[i+1])*(s[i+1] - s[i])
        ]

    Rows = different phases.
    Columns = positions along the submerged member.
    """

    interval_length = (
        distance[:, 1:] - distance[:, :-1]
    )

    average_force = 0.5 * (
        force_per_metre[:, 1:]
        + force_per_metre[:, :-1]
    )

    return np.sum(
        average_force * interval_length,
        axis=1,
    )


# =====================================================================
# 6. CALCULATE THE WHOLE-JACKET FORCE HISTORY
# =====================================================================

def jacket_forces(
    wave_heading,
    current_heading,
    phases_degrees,
    n_points,
):
    wave_direction, current_velocity, k, sigma = wave_parameters(
        wave_heading,
        current_heading,
    )

    phases = np.deg2rad(
        np.atleast_1d(phases_degrees)
    )

    fractions = np.linspace(0.0, 1.0, n_points)

    # Array indices:
    # member group, phase, horizontal component (x or y)
    drag_groups = np.zeros((4, len(phases), 2))
    inertia_groups = np.zeros_like(drag_groups)

    for group, diameter, start, end in MEMBERS:
        start = np.array(start, dtype=float)
        end = np.array(end, dtype=float)

        change = end - start

        length = np.sqrt(np.sum(change**2))
        span = np.sqrt(np.sum(change[:2]**2))

        # Local phase at the start of the member.
        start_phase = (
            phases
            - k * np.sum(start[:2] * wave_direction)
        )

        # Phase difference between the two ends.
        phase_change = (
            k * np.sum(change[:2] * wave_direction)
        )

        wet = submerged_fraction(
            start[2],
            change[2],
            start_phase,
            phase_change,
            fractions,
        )

        if not np.any(wet):
            continue

        # One row per phase; one column per spatial point.
        fraction_along_member = (
            wet[:, None] * fractions
        )

        distance = length * fraction_along_member

        z = (
            start[2]
            + change[2] * fraction_along_member
        )

        theta = (
            start_phase[:, None]
            - phase_change * fraction_along_member
        )

        # -------------------------------------------------------------
        # Wave velocity amplitudes
        # -------------------------------------------------------------

        scale = (
            AMPLITUDE * sigma / np.sinh(k * DEPTH)
        )

        U = scale * np.cosh(k * (z + DEPTH))
        W = scale * np.sinh(k * (z + DEPTH))

        wave_u = U * np.sin(theta)
        w = W * np.cos(theta)

        # Analytical first-order material acceleration.
        wave_acceleration = sigma * U * np.cos(theta)
        az = -sigma * W * np.sin(theta)

        # Add wave and current velocities as vectors.
        ux = (
            current_velocity[0]
            + wave_direction[0] * wave_u
        )

        uy = (
            current_velocity[1]
            + wave_direction[1] * wave_u
        )

        ax = wave_direction[0] * wave_acceleration
        ay = wave_direction[1] * wave_acceleration

        drag_factor = 0.5 * RHO * CD * diameter

        inertia_factor = (
            RHO * CM * np.pi * diameter**2 / 4.0
        )

        # -------------------------------------------------------------
        # Vertical columns
        # -------------------------------------------------------------

        if group == 0:
            speed = np.sqrt(ux**2 + uy**2)

            components = [
                (0, ux, ax),
                (1, uy, ay),
            ]

            for component, velocity, acceleration in components:
                drag_per_metre = (
                    drag_factor * velocity * speed
                )

                inertia_per_metre = (
                    inertia_factor * acceleration
                )

                drag_groups[group, :, component] += integrate_force(
                    drag_per_metre,
                    distance,
                )

                inertia_groups[group, :, component] += integrate_force(
                    inertia_per_metre,
                    distance,
                )

            continue

        # -------------------------------------------------------------
        # Beams and diagonals
        # -------------------------------------------------------------

        along = change[:2] / span

        # First normal direction: horizontal across the member.
        normals = [
            np.array([-along[1], along[0], 0.0])
        ]

        # Second normal direction for diagonals.
        if group in (2, 3):
            sin_alpha = change[2] / length
            cos_alpha = span / length

            normals.append(
                np.array([
                    sin_alpha * along[0],
                    sin_alpha * along[1],
                    -cos_alpha,
                ])
            )

        # The other normal direction of a horizontal beam is vertical.
        # It gives no horizontal force under the component model.
        for nx, ny, nz in normals:
            un = nx * ux + ny * uy + nz * w
            an = nx * ax + ny * ay + nz * az

            drag_per_metre = (
                drag_factor * un * np.abs(un)
            )

            inertia_per_metre = (
                inertia_factor * an
            )

            normal_drag = integrate_force(
                drag_per_metre,
                distance,
            )

            normal_inertia = integrate_force(
                inertia_per_metre,
                distance,
            )

            horizontal_projection = np.array([nx, ny])

            drag_groups[group] += (
                normal_drag[:, None] * horizontal_projection
            )

            inertia_groups[group] += (
                normal_inertia[:, None] * horizontal_projection
            )

    # Add signed forces at the SAME instant.
    drag = np.sum(drag_groups, axis=0)
    inertia = np.sum(inertia_groups, axis=0)

    total = drag + inertia

    resultant = np.sqrt(
        total[:, 0]**2 + total[:, 1]**2
    )

    return {
        "drag": drag,
        "inertia": inertia,
        "total": total,
        "resultant": resultant,
        "groups": drag_groups + inertia_groups,
        "k": k,
        "sigma": sigma,
    }


# =====================================================================
# 7. FIND THE MAXIMUM USING np.argmax
# =====================================================================

def peak_for_pair(
    wave_heading,
    current_heading,
    phases,
    n_points,
):
    phases = np.atleast_1d(phases)

    result = jacket_forces(
        wave_heading,
        current_heading,
        phases,
        n_points,
    )

    index = np.argmax(result["resultant"])

    return {
        "force": result["resultant"][index],
        "wave": wave_heading % 360.0,
        "current": current_heading % 360.0,
        "phase": phases[index] % 360.0,
    }


def scan_grid(
    wave_headings,
    current_headings,
    phases,
    n_points,
):
    peak_map = np.zeros((
        len(wave_headings),
        len(current_headings),
    ))

    phase_map = np.zeros_like(peak_map)

    best = None

    for i, wave in enumerate(wave_headings):
        for j, flow in enumerate(current_headings):
            peak = peak_for_pair(
                wave,
                flow,
                phases,
                n_points,
            )

            peak_map[i, j] = peak["force"]
            phase_map[i, j] = peak["phase"]

            if best is None or peak["force"] > best["force"]:
                best = peak

    return best, peak_map, phase_map


def refine_peak(
    peak,
    heading_width,
    heading_count,
    phase_width,
    phase_count,
    n_points,
):
    """
    Repeat the search on a smaller, finer grid around a candidate peak.
    """

    waves = peak["wave"] + np.linspace(
        -heading_width,
        heading_width,
        heading_count,
    )

    currents = peak["current"] + np.linspace(
        -heading_width,
        heading_width,
        heading_count,
    )

    phases = peak["phase"] + np.linspace(
        -phase_width,
        phase_width,
        phase_count,
    )

    best, peak_map, phase_map = scan_grid(
        waves,
        currents,
        phases,
        n_points,
    )

    return best


def find_maximum():
    # Initial sampling over the full range.
    headings = np.linspace(0.0, 360.0, 31)[:-1]
    phases = np.linspace(0.0, 360.0, 91)[:-1]

    print("Scanning all wave/current heading pairs...")

    coarse_best, peak_map, phase_map = scan_grid(
        headings,
        headings,
        phases,
        65,
    )

    # Find local peaks of the heading map.
    # % wraps around because 360 degrees equals 0 degrees.
    candidates = []

    for i in range(len(headings)):
        for j in range(len(headings)):
            neighbours = []

            for di in [-1, 0, 1]:
                for dj in [-1, 0, 1]:
                    row = (i + di) % len(headings)
                    column = (j + dj) % len(headings)

                    neighbours.append(
                        peak_map[row, column]
                    )

            if peak_map[i, j] >= max(neighbours):
                candidates.append({
                    "wave": headings[i],
                    "current": headings[j],
                    "phase": phase_map[i, j],
                    "force": peak_map[i, j],
                })

    # Refine every candidate instead of assuming the first is best.
    refined = []

    for i, peak in enumerate(candidates):
        print(
            f"Refining candidate {i + 1}/{len(candidates)}..."
        )

        peak = refine_peak(
            peak,
            heading_width=12.0,
            heading_count=13,
            phase_width=12.0,
            phase_count=25,
            n_points=129,
        )

        peak = refine_peak(
            peak,
            heading_width=2.0,
            heading_count=17,
            phase_width=3.0,
            phase_count=25,
            n_points=257,
        )

        peak = refine_peak(
            peak,
            heading_width=0.25,
            heading_count=11,
            phase_width=0.5,
            phase_count=21,
            n_points=513,
        )

        # Final finer phase sampling at this direction pair.
        phases = peak["phase"] + np.linspace(
            -0.05, 0.05, 21
        )

        peak = peak_for_pair(
            peak["wave"],
            peak["current"],
            phases,
            2049,
        )

        refined.append(peak)

    values = np.array([
        peak["force"] for peak in refined
    ])

    best = refined[np.argmax(values)]

    return best, headings, peak_map


# =====================================================================
# 8. RESULTS AND CONVERGENCE CHECK
# =====================================================================

def print_results(best):
    result = jacket_forces(
        best["wave"],
        best["current"],
        best["phase"],
        2049,
    )

    print("\nPART 2(a): LARGEST HORIZONTAL RESULTANT FOUND")

    print(
        f"Force:           {result['resultant'][0] / 1e6:.4f} MN"
    )
    print(f"Wave heading:    {best['wave']:.3f} degrees")
    print(f"Current heading: {best['current']:.3f} degrees")
    print(f"Phase at origin: {best['phase']:.3f} degrees")

    print(
        f"Time:            {best['phase'] / 360.0 * PERIOD:.4f} s"
    )

    print("\nContributions at the SAME INSTANT:")

    print(
        f"{'Contribution':<24}"
        f"{'Fx (MN)':>12}"
        f"{'Fy (MN)':>12}"
    )

    for label in ["drag", "inertia", "total"]:
        fx, fy = result[label][0] / 1e6

        print(
            f"{label.capitalize():<24}"
            f"{fx:>12.4f}"
            f"{fy:>12.4f}"
        )

    for group, name in enumerate(GROUPS):
        fx, fy = result["groups"][group, 0] / 1e6

        print(
            f"{name:<24}"
            f"{fx:>12.4f}"
            f"{fy:>12.4f}"
        )

    print(
        "\nSpatial refinement at these same headings and phase:"
    )
    print("Points/member    Horizontal force (MN)")

    for n_points in [257, 513, 1025, 2049]:
        check = jacket_forces(
            best["wave"],
            best["current"],
            best["phase"],
            n_points,
        )

        print(
            f"{n_points:<17}"
            f"{check['resultant'][0] / 1e6:.7f}"
        )


# =====================================================================
# 9. BASIC PLOTS
# =====================================================================

def plot_results(best, headings, peak_map):
    plt.figure(figsize=(12, 5))

    # Heading comparison
    plt.subplot(1, 2, 1)

    waves, currents = np.meshgrid(
        headings,
        headings,
        indexing="ij",
    )

    plt.pcolormesh(
        currents,
        waves,
        peak_map / 1e6,
        shading="auto",
    )

    plt.colorbar(label="Cycle maximum (MN)")

    plt.scatter(
        best["current"],
        best["wave"],
        color="red",
        edgecolors="white",
    )

    plt.xlabel("Current heading (degrees)")
    plt.ylabel("Wave heading (degrees)")
    plt.title("Independent wave/current directions")

    # Time histories
    plt.subplot(1, 2, 2)

    phases = np.linspace(0.0, 360.0, 721)
    time = phases / 360.0 * PERIOD

    cases = [
        (0.0, 0.0, "Both along +x"),
        (0.0, 90.0, "Wave +x, current +y"),
        (
            best["wave"],
            best["current"],
            "Largest peak found",
        ),
    ]

    for wave, flow, label in cases:
        result = jacket_forces(
            wave,
            flow,
            phases,
            513,
        )

        plt.plot(
            time,
            result["resultant"] / 1e6,
            label=label,
        )

    plt.scatter(
        best["phase"] / 360.0 * PERIOD,
        best["force"] / 1e6,
        color="red",
        zorder=3,
    )

    plt.xlabel("Time (s)")
    plt.ylabel("Horizontal resultant (MN)")
    plt.title("One wave cycle")
    plt.grid(alpha=0.25)
    plt.legend(fontsize=8)

    plt.tight_layout()
    plt.show()


# =====================================================================
# 10. RUN
# =====================================================================

if __name__ == "__main__":
    best, headings, peak_map = find_maximum()

    print_results(best)

    plot_results(
        best,
        headings,
        peak_map,
    )
    
#%%
# Keep this plotting cell import-safe for the viewer exporter.
if __name__ == "__main__":
    # ============================================================
    # 3D LOAD SURFACE
    #
    # X-axis: current direction
    # Y-axis: wave direction
    # Height: maximum horizontal load over a wave cycle
    #
    # Existing results:
    # headings = sampled directions, from 0 up to but excluding 360
    # peak_map[i, j] = maximum load for wave i and current j, in N
    # best = refined maximum and its directions and phase
    # ============================================================

    # 1. Convert force from N to MN
    load_MN = peak_map / 1_000_000

    # 2. Close the surface at 360 degrees.
    #    A heading of 360 degrees is identical to 0 degrees.
    angles = np.append(headings, 360.0)

    n = len(headings)
    surface_load = np.zeros((n + 1, n + 1))

    surface_load[:-1, :-1] = load_MN
    surface_load[-1, :-1] = load_MN[0, :]
    surface_load[:-1, -1] = load_MN[:, 0]
    surface_load[-1, -1] = load_MN[0, 0]

    # 3. Create the two direction grids.
    #    Rows correspond to waves; columns correspond to currents.
    wave_grid, current_grid = np.meshgrid(
        angles, angles, indexing="ij"
    )

    # 4. Draw the surface
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(
        111, projection="3d", computed_zorder=False
    )

    surface = ax.plot_surface(
        current_grid,
        wave_grid,
        surface_load,
        cmap="viridis",
        rstride=1,
        cstride=1,
        edgecolor=(0, 0, 0, 0.12),
        linewidth=0.15,
        alpha=0.95,
        zorder=1
    )

    # 5. Mark the refined maximum
    ax.scatter(
        best["current"],
        best["wave"],
        best["force"] / 1_000_000,
        color="red",
        edgecolors="white",
        s=80,
        depthshade=False,
        zorder=10
    )

    # 6. Display its load, directions and phase
    peak_text = (
        f"Largest load found: {best['force'] / 1_000_000:.3f} MN\n"
        f"Wave direction: {best['wave']:.2f} degrees\n"
        f"Current direction: {best['current']:.2f} degrees\n"
        f"Phase at origin: {best['phase']:.2f} degrees"
    )

    ax.text2D(
        0.02, 0.97,
        peak_text,
        transform=ax.transAxes,
        va="top",
        fontsize=10,
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="none")
    )

    # 7. Labels and viewing angle
    ax.set_xlabel("Current direction (degrees)", labelpad=12)
    ax.set_ylabel("Wave direction (degrees)", labelpad=12)
    ax.set_zlabel("Maximum horizontal load (MN)", labelpad=12)

    ax.set_xlim(0, 360)
    ax.set_ylim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_yticks([0, 90, 180, 270, 360])

    ax.view_init(elev=28, azim=-125)
    ax.set_title("Maximum horizontal load over a wave cycle", pad=25)

    fig.colorbar(
        surface,
        ax=ax,
        shrink=0.6,
        pad=0.12,
        label="Load (MN)"
    )

    fig.text(
        0.5, 0.035,
        "Each direction pair has its own maximizing phase.",
        ha="center",
        fontsize=10
    )

    plt.show()
