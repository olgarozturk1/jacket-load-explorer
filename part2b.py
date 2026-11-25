"""
PART 2(b): VERTICAL LOADING ON THE WHOLE JACKET

Same numerical approach as Part 2(a): NumPy sampling, interpolation,
trapezoidal integration, and a coarse search followed by finer grids.

Assumptions retained:
  - Independent wave and current headings; linear wave-current theory.
  - First-order material acceleration; separate normal-component drag.
  - Drag AND inertia integrated to the instantaneous local wave surface.
  - Linear kinematics continued above SWL; member-centreline wetting.
  - Constant SWL buoyancy for idealised sealed member volumes.
  - Hydrostatic uplift is assumed available, including at column bases.
  - Member-volume overlaps at joints are not subtracted.
  - No axial loading, lift, shielding, diffraction, or impact loading.

Headings: direction of travel, from +x towards +y, in degrees.
Phase: omega*t at the origin. z=0 is SWL; +z is upwards.
Force calculations use N; printed results use MN.
"""

import numpy as np
import matplotlib.pyplot as plt


# 1. INPUTS ------------------------------------------------------------

DEPTH, AMPLITUDE, PERIOD = 40.0, 6.0, 10.0
CURRENT_SPEED = 1.6
RHO, GRAVITY, CD, CM = 1025.0, 9.81, 1.2, 0.8
MASS = 8_000_000.0  # 8000 tonnes, expressed in kg
OMEGA = 2.0 * np.pi / PERIOD
K_SAMPLES = np.linspace(0.0, 0.12, 12001)


# 2. GEOMETRY AND CONSTANT LOADS ---------------------------------------
# Each row: group, diameter (m), start (x,y,z), end (x,y,z).
# Sloping/vertical members start at their lower end.

GROUPS = ["Columns", "Horizontal beams", "x-z diagonals", "y-z diagonals"]

MEMBERS = [
    (0, 4.0, (0, 0, -40), (0, 0, 12)),
    (0, 4.0, (0, 30, -40), (0, 30, 12)),
    (0, 4.0, (30, 0, -40), (30, 0, 12)),
    (0, 4.0, (30, 30, -40), (30, 30, 12)),

    (1, 2.0, (0, 0, -35), (30, 0, -35)),
    (1, 2.0, (0, 30, -35), (30, 30, -35)),
    (1, 2.0, (0, 0, -35), (0, 30, -35)),
    (1, 2.0, (30, 0, -35), (30, 30, -35)),

    (1, 2.0, (0, 0, 7), (30, 0, 7)),
    (1, 2.0, (0, 30, 7), (30, 30, 7)),
    (1, 2.0, (0, 0, 7), (0, 30, 7)),
    (1, 2.0, (30, 0, 7), (30, 30, 7)),

    (2, 1.0, (0, 0, -35), (30, 0, 7)),
    (2, 1.0, (0, 30, -35), (30, 30, 7)),
    (3, 1.0, (0, 30, -35), (0, 0, 7)),
    (3, 1.0, (30, 30, -35), (30, 0, 7)),
]

WEIGHT = MASS * GRAVITY  # positive magnitude, acting downwards
swl_volume = 0.0

for group, diameter, start, end in MEMBERS:
    change = np.array(end, dtype=float) - np.array(start, dtype=float)
    length = np.sqrt(np.sum(change**2))
    if change[2] == 0:
        wet = float(start[2] < 0)
    else:
        wet = np.clip(-start[2] / change[2], 0.0, 1.0)
    swl_volume += np.pi * diameter**2 / 4.0 * length * wet

BUOYANCY = RHO * GRAVITY * swl_volume  # upwards, constant


# 3. WAVE NUMBER AND INTRINSIC FREQUENCY -------------------------------

def wave_parameters(wave_heading, current_heading):
    wave_angle, current_angle = np.deg2rad([wave_heading, current_heading])
    direction = np.array([np.cos(wave_angle), np.sin(wave_angle)])
    current = CURRENT_SPEED * np.array([
        np.cos(current_angle), np.sin(current_angle)
    ])
    current_along_wave = np.sum(current * direction)

    # omega = sqrt(g*k*tanh(k*d)) + k*current_along_wave
    frequency = (
        np.sqrt(GRAVITY * K_SAMPLES * np.tanh(K_SAMPLES * DEPTH))
        + K_SAMPLES * current_along_wave
    )
    if (np.any(np.diff(frequency) <= 0)
            or not frequency[0] <= OMEGA <= frequency[-1]):
        raise ValueError("Reconsider the k sampling interval for these inputs.")

    k = np.interp(OMEGA, frequency, K_SAMPLES)
    sigma = OMEGA - k * current_along_wave
    return direction, current, k, sigma


# 4. SUBMERGED LENGTH AND TRAPEZOIDAL RULE ------------------------------

def submerged_fraction(start_z, rise, start_phase, phase_change, fractions):
    if rise == 0:
        # With this geometry, lower beams are wet and upper beams dry.
        if start_z < -AMPLITUDE:
            return np.ones_like(start_phase)
        if start_z >= AMPLITUDE:
            return np.zeros_like(start_phase)
        raise ValueError("A beam within the wave envelope needs new wetting logic.")

    # This condition guarantees one surface crossing for a rising member.
    if rise <= AMPLITUDE * abs(phase_change):
        raise ValueError("This member may have multiple water-surface crossings.")

    z = start_z + rise * fractions
    wet = np.zeros_like(start_phase)
    for i in range(len(start_phase)):
        eta = AMPLITUDE * np.sin(start_phase[i] - phase_change * fractions)
        wet[i] = np.interp(0.0, z - eta, fractions)
    return wet


def integrate_force(force_per_metre, distance):
    # Rows = phases; columns = positions along the submerged member.
    ds = distance[:, 1:] - distance[:, :-1]
    average = 0.5 * (force_per_metre[:, 1:] + force_per_metre[:, :-1])
    return np.sum(average * ds, axis=1)


# 5. VERTICAL FORCE HISTORY --------------------------------------------

def jacket_forces(wave_heading, current_heading, phases_degrees, n_points):
    direction, current, k, sigma = wave_parameters(wave_heading, current_heading)
    phases = np.deg2rad(np.atleast_1d(phases_degrees))
    fractions = np.linspace(0.0, 1.0, n_points)
    drag_groups = np.zeros((4, len(phases)))
    inertia_groups = np.zeros_like(drag_groups)

    for group, diameter, start, end in MEMBERS:
        # Columns have no vertical cross-flow force. Axial loads are omitted.
        if group == 0:
            continue

        start = np.array(start, dtype=float)
        change = np.array(end, dtype=float) - start
        length = np.sqrt(np.sum(change**2))
        span = np.sqrt(np.sum(change[:2]**2))
        start_phase = phases - k * np.sum(start[:2] * direction)
        phase_change = k * np.sum(change[:2] * direction)

        wet = submerged_fraction(
            start[2], change[2], start_phase, phase_change, fractions
        )
        if not np.any(wet):
            continue

        position = wet[:, None] * fractions
        distance = length * position
        z = start[2] + change[2] * position
        theta = start_phase[:, None] - phase_change * position

        scale = AMPLITUDE * sigma / np.sinh(k * DEPTH)
        U = scale * np.cosh(k * (z + DEPTH))
        W = scale * np.sinh(k * (z + DEPTH))
        wave_u = U * np.sin(theta)
        w = W * np.cos(theta)
        ux = current[0] + direction[0] * wave_u
        uy = current[1] + direction[1] * wave_u

        # First-order material acceleration, using sigma.
        ax = direction[0] * sigma * U * np.cos(theta)
        ay = direction[1] * sigma * U * np.cos(theta)
        az = -sigma * W * np.sin(theta)

        if group == 1:
            nx, ny, nz = 0.0, 0.0, 1.0  # vertical normal to a horizontal beam
        else:
            along = change[:2] / span
            sin_alpha = change[2] / length
            cos_alpha = span / length
            nx, ny, nz = sin_alpha * along[0], sin_alpha * along[1], -cos_alpha

        # The other normal is horizontal and contributes no vertical force
        # under the agreed separate normal-component drag approximation.
        un = nx * ux + ny * uy + nz * w
        an = nx * ax + ny * ay + nz * az
        normal_drag = 0.5 * RHO * CD * diameter * un * np.abs(un)
        normal_inertia = RHO * CM * np.pi * diameter**2 / 4.0 * an

        drag_groups[group] += nz * integrate_force(normal_drag, distance)
        inertia_groups[group] += nz * integrate_force(normal_inertia, distance)

    drag = np.sum(drag_groups, axis=0)
    inertia = np.sum(inertia_groups, axis=0)

    # D and I are upward-positive. C is downward-positive foundation load.
    compression = WEIGHT - BUOYANCY - drag - inertia
    return {
        "drag": drag, "inertia": inertia, "compression": compression,
        "groups": drag_groups + inertia_groups, "k": k, "sigma": sigma,
    }


# 6. SEARCH DIRECTIONS AND PHASE, THEN REFINE ---------------------------

def scan_grid(waves, currents, phases, n_points):
    peak_map = np.zeros((len(waves), len(currents)))
    phase_map = np.zeros_like(peak_map)
    best = None
    for i, wave in enumerate(waves):
        for j, current in enumerate(currents):
            result = jacket_forces(wave, current, phases, n_points)
            index = np.argmax(result["compression"])
            peak_map[i, j] = result["compression"][index]
            phase_map[i, j] = phases[index] % 360.0
            if best is None or peak_map[i, j] > best["force"]:
                best = {"force": peak_map[i, j], "wave": wave % 360.0,
                        "current": current % 360.0, "phase": phase_map[i, j]}
    return best, peak_map, phase_map


def find_maximum():
    headings = np.linspace(0.0, 360.0, 31)[:-1]  # 12-degree spacing
    phases = np.linspace(0.0, 360.0, 91)[:-1]    # 4-degree spacing
    print("Scanning independent wave/current directions and phase...", flush=True)
    _, peak_map, phase_map = scan_grid(headings, headings, phases, 65)

    # Keep every local peak in the heading map; wrap around at 360 degrees.
    candidates = []
    count = len(headings)
    for i in range(count):
        for j in range(count):
            neighbours = [peak_map[(i + di) % count, (j + dj) % count]
                          for di in [-1, 0, 1] for dj in [-1, 0, 1]]
            if peak_map[i, j] >= max(neighbours):
                candidates.append({"wave": headings[i], "current": headings[j],
                                   "phase": phase_map[i, j], "force": peak_map[i, j]})

    refined = []
    for i, peak in enumerate(candidates):
        print(f"Refining candidate {i + 1}/{len(candidates)}...", flush=True)
        # Heading half-width/count, phase half-width/count, spatial points.
        stages = [(12.0, 13, 12.0, 25, 129),
                  (2.0, 17, 3.0, 25, 257),
                  (0.25, 11, 0.5, 21, 513)]
        for hw, hn, pw, pn, points in stages:
            waves = peak["wave"] + np.linspace(-hw, hw, hn)
            currents = peak["current"] + np.linspace(-hw, hw, hn)
            phases = peak["phase"] + np.linspace(-pw, pw, pn)
            peak, _, _ = scan_grid(waves, currents, phases, points)

        phases = peak["phase"] + np.linspace(-0.05, 0.05, 21)
        peak, _, _ = scan_grid([peak["wave"]], [peak["current"]], phases, 2049)
        refined.append(peak)

    best = refined[np.argmax([peak["force"] for peak in refined])]
    return best, headings, peak_map


# 7. CLEAR RESULTS FOR PART 2(b)(i), (ii), AND (iii) --------------------

def print_results(best, phases, cycle):
    result = jacket_forces(best["wave"], best["current"], best["phase"], 2049)
    print("\nPART 2(b): LARGEST DOWNWARD FOUNDATION LOAD FOUND")
    print(f"Compression:      {result['compression'][0] / 1e6:.4f} MN")
    print(f"Wave heading:     {best['wave']:.3f} degrees")
    print(f"Current heading:  {best['current']:.3f} degrees")
    print(f"Phase at origin:  {best['phase']:.3f} degrees")
    print(f"Time:             {best['phase'] / 360.0 * PERIOD:.4f} s")
    print(f"k = {result['k']:.6f} 1/m; sigma = {result['sigma']:.6f} rad/s")

    print("\nCONTRIBUTIONS AT THAT SAME INSTANT (upwards positive)")
    loads = [("Weight", -WEIGHT), ("SWL buoyancy", BUOYANCY),
             ("Vertical drag", result["drag"][0]),
             ("Vertical inertia", result["inertia"][0])]
    for name, value in loads:
        print(f"{name:<22}{value / 1e6:>12.4f} MN")
    print("Downward foundation load = W - B - vertical drag - vertical inertia.")

    print("\nWHOLE-CYCLE STATISTICS AT THESE FIXED HEADINGS (MN)")
    print("Weight, buoyancy, drag and inertia: upwards positive.")
    print("Compression: downwards positive.")
    print(f"{'Contribution':<20}{'Mean':>12}{'Semi-range':>14}{'Minimum':>12}{'Maximum':>12}")
    histories = [("Weight", np.full_like(phases, -WEIGHT)),
                 ("SWL buoyancy", np.full_like(phases, BUOYANCY)),
                 ("Vertical drag", cycle["drag"]),
                 ("Vertical inertia", cycle["inertia"]),
                 ("Compression", cycle["compression"])]
    for name, values in histories:
        # Equally spaced phases exclude the repeated 360-degree endpoint.
        mean = np.mean(values) / 1e6
        half_range = (np.max(values) - np.min(values)) / 2e6
        print(f"{name:<20}{mean:>12.4f}{half_range:>14.4f}"
              f"{np.min(values) / 1e6:>12.4f}{np.max(values) / 1e6:>12.4f}")
    print("Semi-range = (maximum - minimum)/2; the cycle may be asymmetric.")

    print("\nPHASES OF BOTH SIGNED EXTREMES (degrees at the origin)")
    for name, values in histories[2:]:
        print(f"{name:<20}minimum at {phases[np.argmin(values)]:7.2f}; "
              f"maximum at {phases[np.argmax(values)]:7.2f}")
    print("These cycle statistics use 0.25-degree phase spacing.")
    print("The main maximum above uses the finer search.")

    print("\nDYNAMIC VERTICAL LOAD BY MEMBER GROUP AT THE MAIN MAXIMUM")
    print("Drag + inertia only; upwards positive (MN).")
    for group, name in enumerate(GROUPS):
        print(f"{name:<22}{result['groups'][group, 0] / 1e6:>12.4f}")

    print("\nSPATIAL ACCURACY CHECK AT THE SAME HEADINGS AND PHASE")
    for points in [257, 513, 1025, 2049]:
        check = jacket_forces(best["wave"], best["current"], best["phase"], points)
        print(f"{points:4d} points/member: {check['compression'][0] / 1e6:.7f} MN")


# 8. PLOTS -------------------------------------------------------------

def plot_results(best, headings, peak_map, phases, cycle):
    # Repeat the first point only for plotting a closed wave cycle.
    time = np.append(phases, 360.0) / 360.0 * PERIOD
    fig, axes = plt.subplots(3, 1, figsize=(10, 11), layout="constrained")
    for key, label in [("drag", "Vertical drag"), ("inertia", "Vertical inertia")]:
        values = np.append(cycle[key], cycle[key][0]) / 1e6
        axes[0].plot(time, values, label=label)
    axes[0].axhline(0.0, color="black", linewidth=0.7)
    axes[0].set(title="Dynamic forces at the selected headings",
                ylabel="Force (MN), upwards positive", xlabel="Time (s)")

    compression = np.append(cycle["compression"], cycle["compression"][0]) / 1e6
    axes[1].plot(time, compression, label="Foundation compression")
    axes[1].axhline((WEIGHT - BUOYANCY) / 1e6, linestyle="--", color="grey",
                    label="Weight - constant buoyancy")
    axes[1].scatter(best["phase"] / 360.0 * PERIOD, best["force"] / 1e6,
                    color="red", zorder=3, label="Largest compression found")
    axes[1].set(ylabel="Compression (MN), downwards positive", xlabel="Time (s)")
    for ax in axes[:2]:
        ax.grid(alpha=0.25)
        ax.legend()

    surface = axes[2].pcolormesh(headings, headings, peak_map / 1e6, shading="auto")
    axes[2].scatter(best["current"], best["wave"], color="red", edgecolors="white")
    axes[2].set(xlabel="Current heading (degrees)", ylabel="Wave heading (degrees)",
                title="Cycle maximum compression for each direction pair")
    fig.colorbar(surface, ax=axes[2], label="Maximum compression (MN)")
    plt.show()


# 9. RUN ---------------------------------------------------------------

if __name__ == "__main__":
    best, headings, peak_map = find_maximum()
    phases = np.linspace(0.0, 360.0, 1441)[:-1]
    cycle = jacket_forces(best["wave"], best["current"], phases, 513)
    print_results(best, phases, cycle)
    plot_results(best, headings, peak_map, phases, cycle)
