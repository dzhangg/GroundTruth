# GroundTruth

## Question

Does earthquake frequency fall off exponentially with magnitude? If so, at what rate? The Gutenberg-Richter law predicts `log10(N(≥M)) = a − b·M`, where `N(≥M)` is the count of earthquakes at or above magnitude `M`. The slope, the **b-value**, describes how many small quakes accompany each large one (b ≈ 1.0 is typical; a lower b-value means large events make up a bigger share of total seismicity, as in locked subduction zones). This project pulls live global earthquake data and fits that relationship directly from observation.

## Data source

The [USGS Earthquake Hazards Program](https://earthquake.usgs.gov/fdsnws/event/1/) FDSN event API (no API key required). By default the script fetches every magnitude 4.5+ earthquake worldwide over the last 365 days as GeoJSON, which `earthquake_explorer.py` parses into a pandas DataFrame.

## Method

For each magnitude threshold `M` in the observed range (step 0.1), count the earthquakes with magnitude ≥ `M`, giving the cumulative distribution `N(≥M)`. Fit a straight line to `log10(N(≥M))` versus `M` by least squares (`np.polyfit`, degree 1); the fitted slope is `−b`, the intercept is `a`. See `plot_gutenberg_richter()` in [earthquake_explorer.py](earthquake_explorer.py).

## Figure

![Gutenberg-Richter plot](gutenberg_richter.png)

Observed cumulative counts (blue) with the least-squares GR fit (red); the legend reports the fitted b-value for the current data window.

## Running it

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\Activate.ps1
pip install requests pandas numpy matplotlib
python earthquake_explorer.py
```

The script also produces an epicenter map, a depth histogram, a two-region b-value comparison, an aftershock time series, a subduction-zone depth-gradient (Benioff zone) plot, and an interactive Leaflet map (live at https://github.com/dzhangg/GroundTruth/blob/main/earthquakes.geojson). Each run also archives a dated copy of its outputs under `snapshots/` for tracking change over time. Configuration (magnitude cutoff, time window, region filter, transect) lives in the `CONFIG` block at the top of `earthquake_explorer.py`.

## License

Released under the [MIT License](LICENSE).
