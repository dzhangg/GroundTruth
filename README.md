# GroundTruth

Exploring global seismicity with Python. GroundTruth pulls the live USGS
earthquake catalog and turns it into maps, charts, and an interactive web map
that reveal where earthquakes happen, how deep they are, and how often they occur.

## What it does

The script fetches every magnitude 4.5+ earthquake worldwide over a chosen
time window from the USGS FDSN event API (no API key required) and produces:

- **Epicenter map** — every quake plotted on a world map, colored by depth and
  sized by magnitude. Shallow quakes (yellow) trace spreading ridges; deep ones
  (purple) mark subduction zones.
- **Gutenberg-Richter plot** — the magnitude-frequency relationship with an
  estimated b-value from a least-squares fit.
- **Depth histogram** — separates shallow crustal quakes from intermediate and
  deep-focus events at subduction zones.
- **b-value comparison** — overlays GR fits for two configurable regions on one
  plot so you can see how their seismicity distributions differ.
- **Interactive web map** (`index.html`) — Leaflet map with OSM tiles where
  every dot is colored by magnitude, sized by energy release, and opens a popup
  with depth, date, coordinates, and a direct link to the USGS event page.

## Interactive map

The `index.html` file is a self-contained Leaflet map that reads
`earthquakes.geojson`. To use it:

**Locally:**
```bash
python -m http.server
# open http://localhost:8000
On GitHub Pages:

Push this repo to GitHub
Go to Settings → Pages → Source → Deploy from branch → main / (root)
Your map will be live at https://<your-username>.github.io/GroundTruth/
Color key
Color	Magnitude	Class
Green	M < 5.0	Minor
Lime	M 5.0–5.5	Light
Yellow	M 5.5–6.0	Moderate
Orange	M 6.0–6.5	Strong
Red	M 6.5–7.0	Major
Dark red	M ≥ 7.0	Great
Running it

python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\Activate.ps1
pip install requests pandas numpy matplotlib
python earthquake_explorer.py
Output files are saved in the project folder. On the first run, Natural Earth
country boundaries (~14 MB) are downloaded once and cached as countries.geojson.

Configuration
All options live at the top of earthquake_explorer.py:

Variable	Default	Description
MIN_MAGNITUDE	4.5	Minimum magnitude to fetch
DAYS_BACK	365	Days of history to pull
REGION_NAME	None	Short label for a region filter (e.g. "Japan")
MIN_LAT / MAX_LAT	None	Latitude bounds for the region filter
MIN_LON / MAX_LON	None	Longitude bounds for the region filter
COMPARE_REGIONS	Japan vs. South America	Two regions to compare b-values; set to None to skip
Example — zoom into California:


REGION_NAME = "California"
MIN_LAT, MAX_LAT =  32.0,  42.0
MIN_LON, MAX_LON = -124.5, -114.0
Example — compare Japan vs. Iceland:


COMPARE_REGIONS = [
    {"name": "Japan",   "min_lat": 30,  "max_lat": 46, "min_lon": 129, "max_lon": 146},
    {"name": "Iceland", "min_lat": 63,  "max_lat": 67, "min_lon": -25, "max_lon": -13},
]
Built with
Python · requests · pandas · NumPy · matplotlib · Leaflet.js

Data from the USGS Earthquake Hazards Program.

Roadmap
 Region filtering to zoom into a single plate boundary
 Compare b-values between two regions
 Quantify the depth gradient across a subduction zone (Benioff zone)
 Aftershock time-series view (Omori's law)
 Date-stamped output files for tracking change over time


