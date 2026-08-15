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

## Other analyses

The same live dataset feeds a few other views. Each is one function in `earthquake_explorer.py` and one output file; examples below are from a live run.

### Epicenter map

![Epicenter map](map_epicenters.png)

Every fetched event plotted on a world map, colored by depth (bright = shallow, dark = deep) and sized by magnitude. Shallow events trace spreading ridges; the darker clusters mark subduction zones. See `plot_world_map()`.

### Depth histogram

![Depth histogram](depth_histogram.png)

Distribution of focal depths, with dashed lines at the 70 km (crustal/intermediate) and 300 km (intermediate/deep) boundaries. Most seismicity is shallow; the small bump past 500 km is deep-focus subduction activity. See `plot_depth_histogram()`.

### b-value comparison

![Gutenberg-Richter comparison](gr_comparison.png)

Overlays the GR fit for two regions (Japan vs. South America by default) to compare seismicity distributions — set by `COMPARE_REGIONS`. See `plot_gr_comparison()`.

### Aftershock sequence

![Aftershock series](aftershock_series.png)

Finds the largest event in the current dataset, collects aftershocks within 200 km / 90 days, and fits an Omori–Utsu decay law `n(t) = K / (c + t)^p` to the daily rate. See `plot_aftershock_series()`.

### Benioff zone depth gradient

![Benioff zone](benioff_zone.png)

Projects events onto a configurable trench-to-backarc transect (`BENIOFF_TRANSECT`, Tohoku/NE Japan by default) and plots depth vs. along-track distance; the linear fit's slope gives the slab's dip angle. See `plot_benioff_zone()`.

### Interactive map

A Leaflet map reading `earthquakes.geojson`, colored by magnitude with a popup per event — live at https://github.com/dzhangg/GroundTruth/blob/main/earthquakes.geojson.

## Running it

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\Activate.ps1
pip install requests pandas numpy matplotlib
python earthquake_explorer.py
```

Every run also archives a dated copy of its outputs under `snapshots/` for tracking change over time. Configuration (magnitude cutoff, time window, region filter, comparison regions, transect) lives in the `CONFIG` block at the top of `earthquake_explorer.py`.

## Roadmap

- [ ] Pair catalog events with InSAR coseismic interferograms: for M6.5+ events, pull pre/post-event Sentinel-1 SAR pairs from the [Alaska Satellite Facility](https://asf.alaska.edu/) archive, generate a coseismic interferogram, and compare the deformation fringe pattern's location and inferred fault geometry against the catalog epicenter and depth.

## License

Released under the [MIT License](LICENSE).
