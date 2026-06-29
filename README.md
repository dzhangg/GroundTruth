# GroundTruth

Exploring global seismicity with Python. GroundTruth pulls the live USGS
earthquake catalog and turns it into maps and charts that reveal where
earthquakes happen, how deep they are, and how often they occur.

## What it does

The script fetches every magnitude 4.5+ earthquake worldwide over a chosen
time window from the USGS FDSN event API (no API key required) and generates
three figures:

- **Epicenter map** — every quake plotted by location, colored by depth and
  sized by magnitude. The tectonic plate boundaries draw themselves.
- **Gutenberg-Richter plot** — the magnitude-frequency relationship, with an
  estimated b-value from a least-squares fit.
- **Depth distribution** — a histogram that separates shallow crustal quakes
  from the deep events at subduction zones.

## Example output

![Global epicenter map](examples/map_epicenters.png)

A single run of M4.5+ events over the last 365 days (~8,500 earthquakes).
Note how the points trace spreading ridges (thin, shallow, yellow) versus
subduction zones (broad bands deepening inland from yellow to purple).

## Running it

```bash
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\Activate.ps1
pip install requests pandas numpy matplotlib
python earthquake_explorer.py
```

The figures are saved as PNGs in the project folder.

## Configuration

Edit the variables at the top of `earthquake_explorer.py`:

- `MIN_MAGNITUDE` — minimum magnitude to fetch (default 4.5)
- `DAYS_BACK` — how many days of history to pull (default 365)

## Built with

Python, requests, pandas, NumPy, matplotlib. Data from the
[USGS Earthquake Hazards Program](https://earthquake.usgs.gov/fdsnws/event/1/).

## Roadmap

- [ ] Region filtering to zoom into a single plate boundary
- [ ] Compare b-values between two regions
- [ ] Quantify the depth gradient across a subduction zone (Benioff zone)
- [ ] Aftershock time-series view (Omori's law)
- [ ] Date-stamped output files for tracking change over time
