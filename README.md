# Jacket Load Explorer

Part 2(a): horizontal force on the whole jacket under waves and current, with an interactive 3D load surface.

**[Open the public interactive viewer](https://jacket-load-explorer.y655375436.chatgpt.site)**

Drag to rotate, scroll to zoom, choose a viewpoint, and play or select a phase. The height is the horizontal resultant force, not physical elevation on the jacket. The red point marks the largest sampled load at the current phase.

**[Download the AVI](https://jacket-load-explorer.y655375436.chatgpt.site/jacket_load_cycle.avi)** · **[Download an offline interactive copy](https://jacket-load-explorer.y655375436.chatgpt.site/jacket_load_interactive.html)**

## Files

| File | Purpose |
| --- | --- |
| `part2a.py` | The supplied NumPy calculation, grid search, refinement, results and plots. |
| `build_viewer.py` | Recalculates the phase-dependent data using `part2a.py` and builds the interactive HTML. |
| `export_avi.py` | Exports the generated cycle as an AVI with a fixed viewpoint. |
| `viewer/index.html` | Viewer page and labels. |
| `viewer/viewer.js` | Rotation, phase playback, slider and viewpoint controls. |
| `viewer/styles.css` | Viewer layout. |
| `viewer/loads.js` | Saved load data used by the current viewer. |

The calculation functions and input values are the supplied code. The additional 3D plotting cell is guarded so the exporter can import those functions without starting the search or opening figures.

## Run the calculation

Use Python 3.10 or later. Install packages in the same environment used by Spyder or your terminal:

```bash
python -m pip install -r requirements.txt
python part2a.py
```

The search prints its progress, the largest resultant found, its wave/current headings and phase, component forces at that same instant, and a spatial refinement check. If the Matplotlib window blocks the script, close the first figures to continue to the separate 3D surface plot.

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
- Signed member forces are summed at the same instant before taking the horizontal resultant.

The interactive animation uses 12-degree heading spacing, 100 phases per cycle and 129 points along each submerged member. It displays instantaneous load; the static peak surface instead maximizes over phase. The animation's sampled maxima do not replace the finer maximum search in `part2a.py`.

The scientific calculation uses NumPy and Matplotlib. Plotly supplies the rotatable browser chart. The generated Plotly bundle retains its own licence notices; see `THIRD_PARTY_NOTICES.md`.
