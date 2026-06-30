 GroundTruth

Exploring global seismicity with Python. GroundTruth pulls the live USGS earthquake catalog and turns it into maps, charts, and an interactive web map that show where earthquakes happen, how deep they are, and how often they occur.

**Live map:** https://dzhangg.github.io/GroundTruth/
**Repo:** https://github.com/dzhangg/GroundTruth

## What it does

The script fetches every magnitude 4.5+ earthquake worldwide over a configurable time window from the USGS FDSN event API (no API key required), parses the GeoJSON response into a pandas DataFrame, and produces:

- **Epicenter map** (`map_epicenters.png`): every quake plotted on a Natural Earth world map, colored by depth (plasma colormap, shallow = bright, deep = dark, capped at 700 km) and sized by magnitude. Shallow events trace spreading ridges; deep events mark subduction zones.
- **Gutenberg-Richter plot** (`gutenberg_richter.png`): log10(cumulative count) versus magnitude, with a b-value estimated from a least-squares fit.
- **Depth histogram** (`depth_histogram.png`): the depth distribution, with dashed lines at the 70 km (crustal / intermediate) and 300 km (intermediate / deep-focus) boundaries.
- **b-value comparison** (`gr_comparison.png`): overlays GR fits for two configurable regions on one plot so you can compare their seismicity distributions. Defaults to Japan versus South America. Set `COMPARE_REGIONS = None` to skip it.
- **Interactive web map** (`index.html`): a Leaflet map with OpenStreetMap tiles that reads `earthquakes.geojson`. Every marker is colored by magnitude and opens a popup with place, magnitude, depth, date, coordinates, and a direct link to the USGS event page. A legend sits in the bottom-right corner.

The exported `earthquakes.geojson` follows the GitHub Simple Style spec (`marker-color`, `marker-size`) and includes `title` and `description` HTML properties, so it also renders in GitHub's native `.geojson` preview, not just the Leaflet page.

## Interactive map

The live map is served via GitHub Pages at https://dzhangg.github.io/GroundTruth/.

To run it locally:

```bash
python -m http.server
# then open http://localhost:8000
```

To deploy your own copy: push the repo, then go to Settings, Pages, Source, Deploy from branch, `main` / `(root)`.

### Color key (web map, by magnitude)

| Color    | Magnitude | Class    |
| -------- | --------- | -------- |
| Green    | M < 5.0   | Minor    |
| Lime     | M 5.0-5.5 | Light    |
| Yellow   | M 5.5-6.0 | Moderate |
| Orange   | M 6.0-6.5 | Strong   |
| Red      | M 6.5-7.0 | Major    |
| Dark red | M >= 7.0  | Great    |

Note: the web map colors markers by magnitude (palette above). The static epicenter PNG colors points by depth instead.

## Running it

Requires Python 3 (any recent 3.x) and the libraries below.

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\Activate.ps1
pip install requests pandas numpy matplotlib
python earthquake_explorer.py
```

Output files are written to the project folder. On the first run, Natural Earth country boundaries (~14 MB) are downloaded once and cached locally as `countries.geojson` so later runs are instant.

When `REGION_NAME` is set, output filenames are prefixed with the region label (for example `california_map_epicenters.png`) so regional runs do not overwrite the global outputs.

## Configuration

All options live in the CONFIG block at the top of `earthquake_explorer.py`:

| Variable            | Default                | Description                                                       |
| ------------------- | ---------------------- | ----------------------------------------------------------------- |
| `MIN_MAGNITUDE`     | `4.5`                  | Minimum magnitude to fetch                                        |
| `DAYS_BACK`         | `365`                  | Days of history to pull                                           |
| `REGION_NAME`       | `None`                 | Short label for a region filter (used in titles and filenames)    |
| `MIN_LAT`/`MAX_LAT` | `None`                 | Latitude bounds for the optional bounding-box filter              |
| `MIN_LON`/`MAX_LON` | `None`                 | Longitude bounds for the optional bounding-box filter             |
| `COMPARE_REGIONS`   | Japan vs South America | Two regions for the b-value comparison plot; set `None` to skip   |

### Example: zoom into California

```python
REGION_NAME = "California"
MIN_LAT, MAX_LAT =  32.0,  42.0
MIN_LON, MAX_LON = -124.5, -114.0
```

### Example: compare Japan versus Iceland

```python
COMPARE_REGIONS = [
    {"name": "Japan",   "min_lat": 30, "max_lat": 46, "min_lon": 129, "max_lon": 146},
    {"name": "Iceland", "min_lat": 63, "max_lat": 67, "min_lon": -25, "max_lon": -13},
]
```

## Built with

Python, requests, pandas, NumPy, matplotlib, Leaflet.js.

Earthquake data from the [USGS Earthquake Hazards Program](https://earthquake.usgs.gov/fdsnws/event/1/). Country boundaries from [Natural Earth](https://www.naturalearthdata.com/) (public domain).

## Roadmap

- [x] Geographic bounding-box filter to zoom into a single plate boundary
- [x] Compare b-values between two regions
- [ ] Aftershock time series with an Omori-Utsu decay law fit, `n(t) = K / (c + t)^p`
- [ ] Depth gradient across a subduction zone (Benioff zone)
- [ ] Date-stamped output filenames for tracking change over time


