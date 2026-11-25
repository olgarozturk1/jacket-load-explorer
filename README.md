# Jacket Load Explorer

Parts 2(a) and 2(b): horizontal and vertical loading on the whole jacket under waves and current. The interactive 3D viewer illustrates the horizontal load from Part 2(a).

## Start here: Main results — Part 2(a)

**[Read the explanation and main results (PDF)](docs/Hydro%20CW%20Part%202%20-%20a.pdf)**

**[Run this file: part2a.py](part2a.py)**

This is the main calculation for Part 2(a). It finds and prints the largest horizontal load, the corresponding wave and current directions, the phase, and the force contributions from each member group. Run it in Spyder, or follow [the running instructions below](#run-the-calculation).

## Main results — Part 2(b)

**[Run this file: part2b.py](part2b.py)**

This is the main calculation for Part 2(b). It uses NumPy sampling, interpolation, trapezoidal integration and progressively finer search grids, following the same numerical approach as Part 2(a).

For the supplied inputs, the largest downward foundation load found is **53.6927 MN**, including weight, constant SWL buoyancy, vertical drag and vertical inertia. The refined grid gives wave heading **353.25°**, current heading **344.55°**, and reference phase **151.335°**. These are numerical search results under the stated assumptions.

The script prints the contributions at that same instant, cycle means and semi-ranges, both signed extremes and their phases, member-group contributions, and a spatial accuracy check. Cycle statistics are evaluated at the selected heading pair. It also plots the vertical force histories and a map of cycle-maximum compression against wave and current headings.

Vertical drag and inertia are positive upwards. Foundation compression is positive downwards: **C = W − B − F_D,z − F_I,z**, where W and B are positive magnitudes. Semi-range means (maximum − minimum)/2; it need not describe a symmetric oscillation about the mean.

## Interactive viewer and downloads

**[Open the public interactive viewer](https://jacket-load-explorer.y655375436.chatgpt.site)**

Drag to rotate, scroll to zoom, choose a viewpoint, and play or select a phase. The height is the horizontal resultant force, not physical elevation on the jacket. The red point marks the largest sampled load at the current phase.

**[Download the AVI](https://jacket-load-explorer.y655375436.chatgpt.site/jacket_load_cycle.avi)** · **[Download an offline interactive copy](https://jacket-load-explorer.y655375436.chatgpt.site/jacket_load_interactive.html)**

## Files

| File | Purpose |
| --- | --- |
| [Hydro CW Part 2 - a.pdf](docs/Hydro%20CW%20Part%202%20-%20a.pdf) | Explanation of the search, assumptions, main results and spatial refinement check. |
| **[part2a.py](part2a.py)** | **Main results — run this file for Part 2(a).** Includes the maximum-load search, results and plots. |
| **[part2b.py](part2b.py)** | **Main results — run this file for Part 2(b).** Vertical loads, constant buoyancy, maximum foundation compression, cycle statistics and plots. |
| `build_viewer.py` | Recalculates the phase-dependent data using `part2a.py` and builds the interactive HTML. |
| `export_avi.py` | Exports the generated cycle as an AVI with a fixed viewpoint. |
| `viewer/index.html` | Viewer page and labels. |
| `viewer/viewer.js` | Rotation, phase playback, slider and viewpoint controls. |
| `viewer/styles.css` | Viewer layout. |
| `viewer/loads.js` | Saved load data used by the current viewer. |

The Part 2(a) calculation functions and input values are the supplied code. Its additional 3D plotting cell is guarded so the exporter can import those functions without starting the search or opening figures. Both calculation files can also be imported without running their searches.

## Run the calculation

Use Python 3.10 or later. Install packages in the same environment used by Spyder or your terminal:

```bash
python -m pip install -r requirements.txt
python part2a.py
```

For Part 2(b), run this instead after installing the packages:

```bash
python part2b.py
```

In Spyder, open the required file and run the whole script. Both calculations use NumPy and Matplotlib; Part 2(b) does not require SciPy. The searches print progress and their refined results. For Part 2(a), close the first figures to continue to the separate 3D surface plot if the Matplotlib window blocks the script.

## Build and share the interactive viewer

```bash
python build_viewer.py
```

Open `viewer/jacket_load_interactive.html` in a browser. This single file includes the data and plotting library, so it can be sent to someone else and viewed without Python or an internet connection. It is generated locally rather than committed as a duplicate of the source.

For local development after building, you can also serve the viewer folder:

```bash
python -m http.server 8000 --directory viewer
```

Then open `http://localhost:8000`. Editing this repository does not automatically change the existing public Site; publication is a separate step.

## Export an AVI

```bash
python -m pip install -r requirements-video.txt
python build_viewer.py
python export_avi.py
```

The output is `generated/jacket_load_cycle.avi`. If you already ran the builder for the current inputs, it need not be run again. AVI is a recorded view; use the HTML viewer for rotation and zoom.

## Inputs and model

- Water depth 40 m; wave amplitude 6 m; period 10 s; uniform current speed 1.6 m/s.
- Wave and current directions vary independently. Positive headings rotate from +x towards +y; phase is omega times time at the origin.
- Linear wave–current dispersion and kinematics, with first-order material acceleration.
- Drag and inertia are integrated along the instantaneous submerged length, using the local wave surface. Linear profiles are continued above still water level.
- Columns use vector horizontal drag. Beams and diagonals retain the independent normal-component drag approximation in the supplied calculation.
- Member-centreline wetting; no axial drag, shielding or diffraction.
- Signed member forces are summed at the same instant. Part 2(a) takes the horizontal resultant; Part 2(b) combines vertical drag and inertia with weight and buoyancy to calculate downward foundation compression.

Part 2(b) additionally assumes an 8000-tonne structural mass and constant SWL buoyancy based on idealised sealed member volumes. The corresponding hydrostatic uplift is assumed available, including at column bases; overlapping volumes at joints are not subtracted. Weight and buoyancy have zero oscillation amplitude. Vertical cross-flow forces on the columns are zero in this model because axial loading is omitted. Both dynamic force terms still use the instantaneous local water surface, rather than fixed SWL or crest integration limits.

The interactive animation uses 12-degree heading spacing, 100 phases per cycle and 129 points along each submerged member. It displays instantaneous load; the static peak surface instead maximizes over phase. The animation's sampled maxima do not replace the finer maximum search in `part2a.py`.

The scientific calculation uses NumPy and Matplotlib. Plotly supplies the rotatable browser chart. The generated Plotly bundle retains its own licence notices; see `THIRD_PARTY_NOTICES.md`.
