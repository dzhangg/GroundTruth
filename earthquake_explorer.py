"""
earthquake_explorer.py
======================
Downloads the last N days of earthquakes from the USGS FDSN event API,
parses them into a pandas DataFrame, and produces three plots:
  1. A world map of epicenters
  2. A Gutenberg-Richter frequency-magnitude plot
  3. A depth histogram

Run:
    python earthquake_explorer.py
"""

import json
import math
import os
import shutil
from datetime import datetime, timedelta, timezone

import matplotlib
matplotlib.use("Agg")  # must precede pyplot import
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MplPolygon
import numpy as np
import pandas as pd
import requests

# =============================================================================
# CONFIG  ← change these values to adjust the query
# =============================================================================
MIN_MAGNITUDE = 4.5    # only fetch quakes at or above this magnitude
DAYS_BACK     = 365    # how many days of history to fetch

# --- Gutenberg-Richter fit method -----------------------------------------------
# "mle" uses the Aki-Utsu maximum-likelihood b-value estimator (recommended —
# unlike the OLS fit it doesn't require binning the data, isn't skewed by
# how the tail is binned, and has a known standard error). "ols" keeps the
# original least-squares fit on the cumulative log-count curve for comparison.
FIT_METHOD = "mle"   # "mle" | "ols"

# --- Magnitude-type filter -------------------------------------------------------
# The FDSN catalog mixes magnitude scales (Mww, Mb, Ml, Md, ...) computed by
# different methods; fitting one b-value across mixed scales biases the
# result. Only events whose magType is in this list are used for the
# Gutenberg-Richter fit — everything else is still fetched and used by the
# other plots (map, depth histogram, aftershock series, Benioff plot).
# Set MAG_TYPES = None to disable the filter and use every event.
MAG_TYPES = ("mww", "mwc", "mwb", "mwr", "mw")   # moment-magnitude family

# --- World map background -----------------------------------------------------
# GeoJSON country boundaries from Natural Earth (110 m resolution, public domain).
# Downloaded once and cached locally so re-runs are instant.
WORLD_MAP_URL   = "https://raw.githubusercontent.com/datasets/geo-countries/master/data/countries.geojson"
WORLD_MAP_CACHE = "countries.geojson"   # written next to the script on first run

# --- Geographic bounding box (TODO 1) ----------------------------------------
# Set REGION_NAME to a short label (used in titles and filenames).
# Set the four lat/lon bounds to restrict the query to one area.
# Leave all five as None to query the whole world (the default).
#
# Example — Japan:
#   REGION_NAME = "Japan"
#   MIN_LAT, MAX_LAT =  30.0,  46.0
#   MIN_LON, MAX_LON = 129.0, 146.0
#
# Example — California:
#   REGION_NAME = "California"
#   MIN_LAT, MAX_LAT =  32.0,  42.0
#   MIN_LON, MAX_LON = -124.5, -114.0
REGION_NAME = None
MIN_LAT = None
MAX_LAT = None
MIN_LON = None
MAX_LON = None

# --- b-value comparison (TODO 2) ---------------------------------------------
# Two regions to compare on a single Gutenberg-Richter plot.
# Each dict needs: name, min_lat, max_lat, min_lon, max_lon.
# Set COMPARE_REGIONS = None to skip this plot entirely.
#
# Good contrasting pairs to try:
#   Japan (subduction, high stress) vs. Iceland (spreading ridge, low stress)
#   Cascadia vs. the Andes
COMPARE_REGIONS = [
    {"name": "Japan",         "min_lat":  30, "max_lat":  46, "min_lon": 129, "max_lon": 146},
    {"name": "South America", "min_lat": -55, "max_lat":  15, "min_lon": -82, "max_lon": -34},
]

# --- Benioff zone transect (TODO 4) -------------------------------------------
# A trench-to-backarc line used to slice through a subduction zone and plot
# depth vs. distance-along, tracing the dipping Wadati-Benioff zone.
# "start" should sit on the trench (shallow) side, "end" on the backarc
# (deeper) side.  half_width_km is the corridor half-width around the line;
# widen it if too few events fall within it.
# Set BENIOFF_TRANSECT = None to skip this plot.
BENIOFF_TRANSECT = {
    "name":          "NE Japan (Tohoku)",
    "start":         (39.5, 144.0),   # trench side
    "end":           (39.0, 138.0),   # backarc side
    "half_width_km": 100.0,
}


# =============================================================================
# SECTION 1 – BUILD THE API URL AND DOWNLOAD DATA
# =============================================================================

def fetch_earthquakes(min_mag: float, days_back: int,
                      min_lat=None, max_lat=None,
                      min_lon=None, max_lon=None) -> dict:
    """
    Ask the USGS FDSN event API for earthquakes and return the raw GeoJSON dict.

    The API returns a GeoJSON FeatureCollection.  Each Feature looks like:
        {
          "type": "Feature",
          "properties": { "mag": 5.1, "place": "...", "time": 1700000000000, ... },
          "geometry":   { "type": "Point", "coordinates": [lon, lat, depth_km] }
        }

    Docs: https://earthquake.usgs.gov/fdsnws/event/1/
    """
    # Calculate ISO-8601 date strings for 'now' and 'N days ago'
    end_time   = datetime.now(timezone.utc)
    start_time = end_time - timedelta(days=days_back)

    # USGS expects dates in "YYYY-MM-DDTHH:MM:SS" format (UTC)
    fmt = "%Y-%m-%dT%H:%M:%S"

    params = {
        "format":      "geojson",
        "starttime":   start_time.strftime(fmt),
        "endtime":     end_time.strftime(fmt),
        "minmagnitude": min_mag,
        "orderby":     "time",        # newest first
    }

    # Add bounding-box params only when the caller supplied them.
    # The USGS API ignores params that aren't present, so None values
    # must be filtered out rather than passed as the string "None".
    bbox = {
        "minlatitude":  min_lat,
        "maxlatitude":  max_lat,
        "minlongitude": min_lon,
        "maxlongitude": max_lon,
    }
    params.update({k: v for k, v in bbox.items() if v is not None})

    base_url = "https://earthquake.usgs.gov/fdsnws/event/1/query"

    region_label = REGION_NAME or "worldwide"
    print(f"[1/4] Fetching M{min_mag}+ earthquakes ({region_label}) from "
          f"{start_time.date()} to {end_time.date()} …", flush=True)

    # requests.get() sends an HTTP GET request; raise_for_status() turns
    # non-200 responses into Python exceptions so we hear about failures.
    response = requests.get(base_url, params=params, timeout=60)
    response.raise_for_status()

    data = response.json()          # parse the JSON body into a Python dict
    count = data["metadata"]["count"]
    print(f"    → received {count} events")
    return data


# =============================================================================
# SECTION 2 – PARSE GEOJSON INTO A PANDAS DATAFRAME
# =============================================================================

def parse_to_dataframe(geojson: dict) -> pd.DataFrame:
    """
    Walk the list of GeoJSON Features and pull out the six columns we care about.

    pandas DataFrames are like spreadsheets in memory: each column is a Series
    of values, and you can do fast math across every row at once.
    """
    print("[2/4] Parsing GeoJSON into DataFrame …", flush=True)

    rows = []
    for feature in geojson["features"]:
        props = feature["properties"]
        coords = feature["geometry"]["coordinates"]  # [lon, lat, depth_km]

        # epoch milliseconds → Python datetime (UTC-aware)
        epoch_ms = props["time"]
        dt = datetime.fromtimestamp(epoch_ms / 1000.0, tz=timezone.utc)

        rows.append({
            "time":      dt,
            "magnitude": props["mag"],
            "mag_type":  props.get("magType") or "unknown",
            "place":     props["place"],
            "longitude": coords[0],
            "latitude":  coords[1],
            "depth_km":  coords[2],
            "usgs_id":   feature.get("id", ""),
            "url":       props.get("url", ""),
        })

    df = pd.DataFrame(rows)

    # Convert types so later math works cleanly
    df["magnitude"] = pd.to_numeric(df["magnitude"], errors="coerce")
    df["depth_km"]  = pd.to_numeric(df["depth_km"],  errors="coerce")

    # Drop rows where essential values are missing
    df = df.dropna(subset=["magnitude", "latitude", "longitude", "depth_km"])

    print(f"    → DataFrame has {len(df)} rows × {len(df.columns)} columns")
    return df


def filter_by_magnitude_type(df: pd.DataFrame, mag_types) -> pd.DataFrame:
    """
    Restrict to events whose magType is in mag_types (case-insensitive).

    Used only for the Gutenberg-Richter fit, since mixing magnitude scales
    (Mww, Mb, Ml, Md, ...) biases the b-value; the API has no server-side
    parameter for this, so it's applied here after parsing. Pass
    mag_types=None to skip filtering and keep every event.
    """
    if mag_types is None:
        return df

    allowed = {t.lower() for t in mag_types}
    keep = df["mag_type"].str.lower().isin(allowed)
    dropped = df.loc[~keep]

    if len(dropped):
        breakdown = ", ".join(
            f"{mtype}={count}" for mtype, count in dropped["mag_type"].value_counts().items()
        )
        print(f"    → dropped {len(dropped)} event(s) outside MAG_TYPES "
              f"{tuple(sorted(allowed))}: {breakdown}")

    return df[keep].copy()


# =============================================================================
# SECTION 3 – PRINT THE LARGEST EVENT
# =============================================================================

def print_largest_event(df: pd.DataFrame) -> None:
    """Find the row with the maximum magnitude and print a one-line summary."""
    print("[3/4] Finding largest event …", flush=True)

    # idxmax() returns the integer index of the row with the largest value
    idx   = df["magnitude"].idxmax()
    event = df.loc[idx]

    # Format the UTC datetime nicely (strip the timezone info for readability)
    date_str = event["time"].strftime("%Y-%m-%d")

    print(f"\n{'='*55}")
    print(f"  LARGEST EVENT in the last {DAYS_BACK} days")
    print(f"  Magnitude : {event['magnitude']}")
    print(f"  Place     : {event['place']}")
    print(f"  Date      : {date_str}")
    print(f"{'='*55}\n")


# =============================================================================
# SECTION 4 – PLOTS
# =============================================================================

def _load_world_geojson() -> dict:
    """
    Return Natural Earth country boundaries as a GeoJSON dict.

    On the first call the file is downloaded from WORLD_MAP_URL and saved to
    WORLD_MAP_CACHE.  Every subsequent call reads from that local file so the
    script works offline and doesn't re-fetch on every run.
    """
    if os.path.exists(WORLD_MAP_CACHE):
        with open(WORLD_MAP_CACHE) as f:
            return json.load(f)

    print("    downloading world boundaries (cached for future runs) …", flush=True)
    resp = requests.get(WORLD_MAP_URL, timeout=30)
    resp.raise_for_status()
    with open(WORLD_MAP_CACHE, "w") as f:
        f.write(resp.text)
    return resp.json()


def _world_patches(geojson: dict) -> list:
    """
    Convert a GeoJSON FeatureCollection of country polygons into a list of
    matplotlib Polygon patches, one patch per polygon ring.

    GeoJSON geometry types we handle:
      Polygon      → coordinates = [ outer_ring, hole1, hole2, ... ]
      MultiPolygon → coordinates = [ [outer, holes…], [outer, holes…], … ]

    We only draw the outer ring of each polygon (holes are small islands cut
    from large landmasses and aren't visible at this resolution).
    """
    patches = []
    for feature in geojson["features"]:
        geom = feature["geometry"]
        # Normalise to a list of polygon-coord-groups regardless of geometry type
        if geom["type"] == "Polygon":
            polygon_list = [geom["coordinates"]]
        else:  # MultiPolygon
            polygon_list = geom["coordinates"]

        for poly in polygon_list:
            # poly[0] is the outer ring; poly[1:] are holes — we skip holes
            exterior = np.array(poly[0])
            patches.append(MplPolygon(exterior, closed=True))
    return patches


def plot_world_map(df: pd.DataFrame, filename: str = "map_epicenters.png") -> None:
    """
    Scatter plot of earthquake epicenters on a lon/lat grid.

    Colour encodes depth (shallow = yellow, deep = purple).
    Point size encodes magnitude (bigger quake → bigger dot).
    """
    fig, ax = plt.subplots(figsize=(14, 7))

    # Zoom to the bounding box if one is configured, otherwise show the full globe.
    # Add 2° of padding so epicenters near the edge aren't clipped.
    pad = 2
    ax.set_xlim(
        (MIN_LON - pad) if MIN_LON is not None else -180,
        (MAX_LON + pad) if MAX_LON is not None else  180,
    )
    ax.set_ylim(
        (MIN_LAT - pad) if MIN_LAT is not None else -90,
        (MAX_LAT + pad) if MAX_LAT is not None else  90,
    )
    ax.set_facecolor("#c8e6f5")   # light blue = ocean

    # Draw land polygons as a single PatchCollection (much faster than one-by-one).
    # facecolor = parchment tan for land; edgecolor = thin grey country borders.
    world_geojson = _load_world_geojson()
    land = PatchCollection(
        _world_patches(world_geojson),
        facecolor="#e8dfc8",
        edgecolor="#999999",
        linewidth=0.3,
        zorder=1,   # behind the earthquake scatter (zorder=2 default for scatter)
    )
    ax.add_collection(land)

    # Scale marker sizes: area proportional to magnitude^3 keeps big quakes visible
    # np.clip prevents very small or negative magnitudes from causing issues
    sizes = (np.clip(df["magnitude"], 4, 10) ** 3) * 0.5

    scatter = ax.scatter(
        df["longitude"],
        df["latitude"],
        c=df["depth_km"],           # colour axis = depth
        s=sizes,                    # marker size
        cmap="plasma_r",            # shallow (small depth) → bright, deep → dark
        alpha=0.5,
        linewidths=0,
        vmin=0,
        vmax=700,                   # most subduction-zone quakes are < 700 km deep
    )

    # Colourbar explains what the colours mean
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Depth (km)", fontsize=11)

    ax.set_xlabel("Longitude", fontsize=11)
    ax.set_ylabel("Latitude", fontsize=11)
    region_label = f" — {REGION_NAME}" if REGION_NAME else ""
    ax.set_title(
        f"M{MIN_MAGNITUDE}+ Earthquake Epicenters{region_label} — Last {DAYS_BACK} Days\n"
        f"({len(df):,} events)  |  Colour = depth, Size ∝ magnitude",
        fontsize=13,
    )
    ax.grid(color="white", linewidth=0.4, alpha=0.6)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"    saved → {filename}")


def plot_gutenberg_richter(
    df: pd.DataFrame, filename: str = "gutenberg_richter.png"
) -> dict:
    """
    Gutenberg-Richter (GR) frequency-magnitude plot.

    The GR relation states:   log10( N(≥M) ) = a − b·M
    where N(≥M) is the number of earthquakes with magnitude ≥ M.  Plots both
    the cumulative distribution and the non-cumulative (per 0.1-wide bin)
    counts, and fits the b-value from the cumulative curve using FIT_METHOD:

      "mle"  Aki-Utsu maximum-likelihood estimator (default), computed
             directly from the event magnitudes at or above the completeness
             threshold Mc = MIN_MAGNITUDE — not from the binned counts —
             with the Utsu (1965) correction for 0.1-wide binning:
                 b_hat = log10(e) / (mean(M) - (Mc - dM/2))
      "ols"  the original least-squares fit of log10(cumulative count)
             vs. M, kept for comparison.

    Either way, uncertainty is reported as the Shi & Bolt (1982) standard
    error:
        sigma_b = 2.30 * b^2 * sqrt( sum((M_i - mean(M))^2) / (n * (n - 1)) )

    Returns a dict of fit statistics (method, mc, b, sigma_b, a, n).
    """
    dM = 0.1                # magnitude binning interval
    # NOTE: Mc reuses the query's magnitude floor, which assumes the MAG_TYPES
    # subset is itself complete down to MIN_MAGNITUDE. It usually isn't — smaller
    # events are less likely to have a moment-magnitude solution at all, so the
    # filtered catalog thins out below its true completeness magnitude, which
    # biases b downward. Watch for a rollover in the non-cumulative curve below
    # the fit line: if you see one, raise Mc (or MIN_MAGNITUDE) until it's gone.
    mc = MIN_MAGNITUDE      # completeness threshold

    # Build a range of magnitude bins from the data minimum to maximum
    m_min = math.floor(df["magnitude"].min() * 10) / 10   # round down to 0.1
    m_max = math.ceil(df["magnitude"].max()  * 10) / 10   # round up   to 0.1
    magnitudes = np.arange(m_min, m_max + dM, dM)          # step = dM

    # Cumulative counts: for each magnitude threshold M, count quakes with mag ≥ M
    cum_counts = np.array(
        [(df["magnitude"] >= m).sum() for m in magnitudes], dtype=float
    )
    valid       = cum_counts > 0   # keep only bins with ≥1 quake so log10 doesn't blow up
    mag_valid   = magnitudes[valid]
    log10_cum   = np.log10(cum_counts[valid])

    # Non-cumulative counts: events falling in each 0.1-wide bin (not fitted —
    # shown only to compare against the cumulative curve used for the fit)
    bin_edges = np.append(magnitudes, magnitudes[-1] + dM)
    incr_counts, _ = np.histogram(df["magnitude"], bins=bin_edges)
    incr_valid = incr_counts > 0
    mag_incr   = magnitudes[incr_valid]
    log10_incr = np.log10(incr_counts[incr_valid])

    # --- b-value fit -------------------------------------------------------
    events_at_mc = df.loc[df["magnitude"] >= mc, "magnitude"]
    n      = len(events_at_mc)
    mean_m = events_at_mc.mean()

    if FIT_METHOD == "mle":
        b_value = math.log10(math.e) / (mean_m - (mc - dM / 2))
        a_value = math.log10(n) + b_value * mc   # anchors the line at (mc, log10 n)
    elif FIT_METHOD == "ols":
        # np.polyfit returns [slope, intercept]; slope is −b in the GR relation.
        coeffs  = np.polyfit(mag_valid, log10_cum, 1)
        b_value = -coeffs[0]
        a_value =  coeffs[1]
    else:
        raise ValueError(f"Unknown FIT_METHOD: {FIT_METHOD!r} (expected 'mle' or 'ols')")

    # Shi & Bolt (1982) standard error of the b-value
    sigma_b = 2.30 * b_value ** 2 * math.sqrt(
        np.sum((events_at_mc - mean_m) ** 2) / (n * (n - 1))
    )

    # Evaluate the fitted line at a smooth set of magnitude values
    fit_x = np.linspace(mag_valid.min(), mag_valid.max(), 200)
    fit_y = a_value - b_value * fit_x

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(mag_valid, log10_cum, color="steelblue", s=30, marker="o",
               label="Cumulative N(≥M)", zorder=3)
    ax.scatter(mag_incr, log10_incr, color="seagreen", s=30, marker="^",
               label="Non-cumulative n(M), 0.1-wide bins", zorder=2)

    ax.plot(fit_x, fit_y, color="tomato", linewidth=2,
            label=f"{FIT_METHOD.upper()} fit  (b = {b_value:.2f} ± {sigma_b:.2f}, n = {n:,})")

    ax.set_xlabel("Magnitude (M)", fontsize=12)
    ax.set_ylabel("log₁₀ Count", fontsize=12)
    region_label = f" — {REGION_NAME}" if REGION_NAME else ""
    ax.set_title(
        f"Gutenberg-Richter Plot{region_label} — Last {DAYS_BACK} Days (M ≥ {MIN_MAGNITUDE})",
        fontsize=13,
    )
    ax.legend(fontsize=10)
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"    saved → {filename}")

    return {
        "method":  FIT_METHOD,
        "mc":      mc,
        "b":       b_value,
        "sigma_b": sigma_b,
        "a":       a_value,
        "n":       n,
    }


def plot_depth_histogram(
    df: pd.DataFrame, filename: str = "depth_histogram.png"
) -> None:
    """
    Histogram of earthquake focal depths.

    Most quakes cluster at shallow depths (< 70 km = crustal),
    with a secondary population of intermediate-depth (70–300 km)
    and deep-focus (> 300 km) events in subduction zones.
    """
    fig, ax = plt.subplots(figsize=(9, 5))

    # np.histogram_bin_edges would also work; here we use fixed 10 km bins
    bins = np.arange(0, df["depth_km"].max() + 20, 20)   # 20 km wide bins

    ax.hist(
        df["depth_km"],
        bins=bins,
        color="mediumseagreen",
        edgecolor="white",
        linewidth=0.4,
    )

    # Mark the three seismological depth zones with vertical dashed lines
    ax.axvline(70,  color="orange", linestyle="--", linewidth=1.2,
               label="Crustal / Intermediate boundary (70 km)")
    ax.axvline(300, color="tomato", linestyle="--", linewidth=1.2,
               label="Intermediate / Deep boundary (300 km)")

    ax.set_xlabel("Depth (km)", fontsize=12)
    ax.set_ylabel("Number of Earthquakes", fontsize=12)
    region_label = f" — {REGION_NAME}" if REGION_NAME else ""
    ax.set_title(
        f"Earthquake Depth Distribution{region_label} — Last {DAYS_BACK} Days (M ≥ {MIN_MAGNITUDE})",
        fontsize=13,
    )
    ax.legend(fontsize=10)
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"    saved → {filename}")


def _transect_project(lat, lon, start: tuple, end: tuple):
    """
    Project point(s) onto the trench-to-backarc line defined by start/end
    (each a (lat, lon) tuple), using a flat-Earth (equirectangular)
    approximation referenced to `start` — accurate enough for the few-hundred
    km scale of a subduction-zone transect.

    lat/lon may be scalars or numpy arrays.  Returns (along_km, perp_km, length_km):
      along_km  distance from `start`, measured along the transect direction
      perp_km   signed perpendicular distance from the transect line
      length_km total length of the transect (start to end)
    """
    lat0, lon0 = start
    lat1, lon1 = end
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * math.cos(math.radians(lat0))

    # Transect direction as a unit vector in local km-space
    ex = (lon1 - lon0) * km_per_deg_lon
    ey = (lat1 - lat0) * km_per_deg_lat
    length = math.hypot(ex, ey)
    ux, uy = ex / length, ey / length

    # Each point's position in the same local km-space
    px = (np.asarray(lon) - lon0) * km_per_deg_lon
    py = (np.asarray(lat) - lat0) * km_per_deg_lat

    along = px * ux + py * uy               # dot product with the unit vector
    perp  = px * uy - py * ux               # cross product = signed offset from the line
    return along, perp, length


def plot_benioff_zone(
    df: pd.DataFrame, transect: dict, filename: str = "benioff_zone.png"
) -> None:
    """
    Depth-gradient plot across a subduction zone (Wadati-Benioff zone).

    Projects every earthquake onto a trench-to-backarc transect (see the
    BENIOFF_TRANSECT config) and plots along-track distance vs. depth for
    events within `half_width_km` of the line.  A straight-line fit to that
    band estimates the slab's dip angle.
    """
    along, perp, length = _transect_project(
        df["latitude"].values, df["longitude"].values,
        transect["start"], transect["end"],
    )
    corridor = df.assign(along_km=along, perp_km=perp)
    corridor = corridor[
        (corridor["perp_km"].abs() <= transect["half_width_km"]) &
        (corridor["along_km"] >= 0) &
        (corridor["along_km"] <= length)
    ]

    name = transect.get("name", "")
    if len(corridor) < 5:
        print(f"    Only {len(corridor)} event(s) within {transect['half_width_km']:.0f} km "
              f"of the {name or 'transect'} line — skipping Benioff zone plot")
        return

    fig, ax = plt.subplots(figsize=(9, 6))

    sc = ax.scatter(
        corridor["along_km"], corridor["depth_km"],
        c=corridor["magnitude"], cmap="viridis",
        s=(np.clip(corridor["magnitude"], 4, 9) ** 2) * 1.5,
        alpha=0.75, linewidths=0,
    )
    cbar = plt.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Magnitude", fontsize=11)

    # Fit a line to depth vs. along-track distance; slope's arctangent is the dip.
    coeffs = np.polyfit(corridor["along_km"], corridor["depth_km"], 1)
    dip_deg = math.degrees(math.atan(coeffs[0]))
    fit_x = np.linspace(corridor["along_km"].min(), corridor["along_km"].max(), 200)
    fit_y = np.polyval(coeffs, fit_x)
    ax.plot(fit_x, fit_y, color="tomato", linewidth=2,
            label=f"Linear fit  (dip ≈ {dip_deg:.0f}°)")

    ax.invert_yaxis()  # depth increases downward on the page
    ax.set_xlabel("Distance from trench-side transect start (km)", fontsize=12)
    ax.set_ylabel("Depth (km)", fontsize=12)
    title_region = f" — {name}" if name else ""
    ax.set_title(
        f"Wadati-Benioff Zone Depth Gradient{title_region}\n"
        f"({len(corridor):,} events within {transect['half_width_km']:.0f} km of transect)",
        fontsize=13,
    )
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"    saved → {filename}")


def _haversine(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance in km between two lat/lon points."""
    R = 6371.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return R * 2 * math.asin(math.sqrt(a))


def plot_aftershock_series(
    df: pd.DataFrame,
    radius_km: float = 200.0,
    days_after: int = 90,
    filename: str = "aftershock_series.png",
) -> None:
    """
    Find the largest event in df, collect aftershocks within radius_km and
    days_after of it, then produce a two-panel figure:

      Top    — magnitude vs. time scatter (colour = depth)
      Bottom — daily aftershock rate bar chart with Omori–Utsu decay fit

    Omori–Utsu law:  n(t) = K / (c + t)^p
      p ≈ 1   typical; higher p → faster decay
      c       small offset that prevents a singularity at t = 0 (fixed at 0.1 day)
      K       productivity constant
    Fit is obtained by log-linear regression on log(n) vs. log(c + t).
    """
    # Identify mainshock
    idx   = df["magnitude"].idxmax()
    ms    = df.loc[idx]
    ms_lat, ms_lon = ms["latitude"], ms["longitude"]
    ms_time = ms["time"]
    ms_mag  = ms["magnitude"]
    ms_place = ms["place"]

    # Filter: events after the mainshock, within radius_km, within days_after days
    cutoff = ms_time + timedelta(days=days_after)
    after  = df[(df["time"] > ms_time) & (df["time"] <= cutoff)].copy()
    after["dist_km"] = after.apply(
        lambda r: _haversine(ms_lat, ms_lon, r["latitude"], r["longitude"]), axis=1
    )
    after = after[after["dist_km"] <= radius_km].copy()

    if len(after) < 5:
        print(
            f"    Only {len(after)} aftershock(s) found within {radius_km} km of "
            f"M{ms_mag:.1f} mainshock — skipping aftershock plot"
        )
        return

    after["days_after"] = (after["time"] - ms_time).dt.total_seconds() / 86400.0

    # Daily rate for Omori–Utsu fit
    day_edges  = np.arange(0, days_after + 1, 1)
    day_ctrs   = day_edges[:-1] + 0.5
    daily_cnt, _ = np.histogram(after["days_after"], bins=day_edges)

    valid = daily_cnt > 0
    c = 0.1  # fixed offset (days)
    log_n = np.log(daily_cnt[valid].astype(float))
    log_t = np.log(day_ctrs[valid] + c)
    coeffs  = np.polyfit(log_t, log_n, 1)
    p_val   = -coeffs[0]
    K_val   = math.exp(coeffs[1])

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 8), gridspec_kw={"height_ratios": [2, 1]}
    )

    # — top: scatter of individual aftershocks —
    sc = ax1.scatter(
        after["days_after"], after["magnitude"],
        c=after["depth_km"], cmap="plasma_r",
        s=(np.clip(after["magnitude"], 4, 9) ** 2) * 2,
        alpha=0.65, linewidths=0, vmin=0, vmax=300,
    )
    cbar = plt.colorbar(sc, ax=ax1, fraction=0.02, pad=0.02)
    cbar.set_label("Depth (km)", fontsize=10)

    ax1.axhline(ms_mag, color="crimson", linestyle="--", linewidth=1.4,
                label=f"Mainshock M{ms_mag:.1f}")
    ax1.set_ylabel("Magnitude", fontsize=11)
    ax1.set_title(
        f"Aftershock Sequence — M{ms_mag:.1f} at {ms_place}\n"
        f"({len(after):,} aftershocks within {radius_km:.0f} km over {days_after} days)",
        fontsize=13,
    )
    ax1.legend(fontsize=10)
    ax1.grid(True, linestyle="--", alpha=0.4)
    ax1.set_xlim(0, days_after)

    # — bottom: daily rate + Omori–Utsu —
    ax2.bar(day_ctrs, daily_cnt, width=0.9, color="steelblue",
            alpha=0.65, label="Daily aftershock count")

    t_fit = np.linspace(0.1, days_after, 500)
    n_fit = K_val / (c + t_fit) ** p_val
    ax2.plot(t_fit, n_fit, color="tomato", linewidth=2,
             label=f"Omori–Utsu fit  (p = {p_val:.2f},  K = {K_val:.1f})")

    ax2.set_xlabel("Days after mainshock", fontsize=11)
    ax2.set_ylabel("Aftershocks / day", fontsize=11)
    ax2.set_xlim(0, days_after)
    ax2.legend(fontsize=10)
    ax2.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"    saved → {filename}")


def plot_gr_comparison(
    df_a: pd.DataFrame, label_a: str,
    df_b: pd.DataFrame, label_b: str,
    filename: str = "gr_comparison.png",
) -> None:
    """
    Overlay the Gutenberg-Richter fits for two regions on one plot.

    The b-value (slope of the log10-count vs. magnitude line) reveals how a
    region's seismicity is distributed across magnitudes:
      b ≈ 1.0  global average
      b > 1.0  more small quakes relative to large ones (e.g. volcanic, spreading)
      b < 1.0  relatively more large quakes (e.g. locked subduction zones)

    A lower b-value doesn't mean a region is "more dangerous" by itself, but it
    does mean large events make up a bigger share of the total seismic activity.
    """
    fig, ax = plt.subplots(figsize=(9, 6))

    palette = [("steelblue", "o"), ("tomato", "s")]   # colour + marker shape

    for (df, label), (color, marker) in zip(
        [(df_a, label_a), (df_b, label_b)], palette
    ):
        m_min = math.floor(df["magnitude"].min() * 10) / 10
        m_max = math.ceil( df["magnitude"].max() * 10) / 10
        magnitudes = np.arange(m_min, m_max + 0.1, 0.1)

        counts = np.array(
            [(df["magnitude"] >= m).sum() for m in magnitudes], dtype=float
        )
        valid        = counts > 0
        mag_valid    = magnitudes[valid]
        log10_count  = np.log10(counts[valid])

        coeffs  = np.polyfit(mag_valid, log10_count, 1)
        b_value = -coeffs[0]

        fit_x = np.linspace(mag_valid.min(), mag_valid.max(), 200)
        fit_y = np.polyval(coeffs, fit_x)

        ax.scatter(mag_valid, log10_count,
                   color=color, marker=marker, s=30, alpha=0.7, zorder=3)
        ax.plot(fit_x, fit_y, color=color, linewidth=2,
                label=f"{label}  —  b = {b_value:.2f}  (n = {len(df):,})")

    ax.set_xlabel("Magnitude (M)", fontsize=12)
    ax.set_ylabel("log₁₀ N(≥M)", fontsize=12)
    ax.set_title(
        f"Gutenberg-Richter Comparison — Last {DAYS_BACK} Days (M ≥ {MIN_MAGNITUDE})\n"
        "Lower b-value → relatively more large quakes",
        fontsize=13,
    )
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"    saved → {filename}")


# =============================================================================
# SECTION 5 – INTERACTIVE GEOJSON (GitHub / Azure Maps viewer)
# =============================================================================

def save_geojson(df: pd.DataFrame, filename: str = "earthquakes.geojson") -> None:
    """
    Write the DataFrame to a GeoJSON FeatureCollection file.

    GitHub automatically renders .geojson files as an interactive map using
    Azure Maps (OSM tiles).  Clicking any dot opens a popup built from the
    'title' and 'description' properties (HTML is supported).

    Colors follow magnitude (intensity), matching the USGS hazard palette:
      M < 5.0   → green    #00b300  minor
      M 5.0–5.5 → lime     #80cc00  light
      M 5.5–6.0 → yellow   #ffcc00  moderate
      M 6.0–6.5 → orange   #ff8000  strong
      M 6.5–7.0 → red      #ff0000  major
      M ≥ 7.0   → dark red #990000  great

    A single "Legend" marker is appended at the end so viewers can decode
    the colour scale without leaving the map.
    """
    def mag_color(mag: float) -> str:
        if mag < 5.0:  return "#00b300"
        if mag < 5.5:  return "#80cc00"
        if mag < 6.0:  return "#ffcc00"
        if mag < 6.5:  return "#ff8000"
        if mag < 7.0:  return "#ff0000"
        return "#990000"

    def mag_symbol_size(mag: float) -> str:
        if mag < 5.5:  return "small"
        if mag < 7.0:  return "medium"
        return "large"

    features = []
    for _, row in df.iterrows():
        mag   = row["magnitude"]
        depth = row["depth_km"]
        lat   = row["latitude"]
        lon   = row["longitude"]
        place = row["place"]
        date  = row["time"].strftime("%Y-%m-%d %H:%M UTC")
        url   = row.get("url", "")

        lat_dir = "N" if lat >= 0 else "S"
        lon_dir = "E" if lon >= 0 else "W"

        # 'description' is rendered as HTML inside the GitHub popup
        desc = (
            f"<strong>M{mag:.1f}</strong> — depth {depth:.0f} km<br>"
            f"{date}<br>"
            f"{abs(lat):.3f}°{lat_dir}, {abs(lon):.3f}°{lon_dir}<br>"
        )
        if url:
            desc += f'<a href="{url}">&#8594; View on USGS</a>'

        features.append({
            "type": "Feature",
            "geometry": {
                "type":        "Point",
                "coordinates": [lon, lat],
            },
            "properties": {
                # 'title' and 'description' used by GitHub's native viewer
                "title":        place,
                "description":  desc,
                # Explicit fields used by index.html (Leaflet)
                "magnitude":    mag,
                "depth_km":     round(depth, 1),
                "date":         row["time"].strftime("%Y-%m-%d %H:%M UTC"),
                "url":          url,
                # GitHub Simple Style — colour and size
                "marker-color": mag_color(mag),
                "marker-size":  mag_symbol_size(mag),
            },
        })

    # Legend point — placed in the South Pacific (open ocean, always on screen
    # for a global view).  Click it to see the full colour key.
    legend_desc = (
        "<strong>Magnitude colour key</strong><br>"
        '<span style="color:#00b300">&#9679;</span> M &lt; 5.0 &nbsp; Minor<br>'
        '<span style="color:#80cc00">&#9679;</span> M 5.0–5.5 &nbsp; Light<br>'
        '<span style="color:#ffcc00">&#9679;</span> M 5.5–6.0 &nbsp; Moderate<br>'
        '<span style="color:#ff8000">&#9679;</span> M 6.0–6.5 &nbsp; Strong<br>'
        '<span style="color:#ff0000">&#9679;</span> M 6.5–7.0 &nbsp; Major<br>'
        '<span style="color:#990000">&#9679;</span> M ≥ 7.0 &nbsp;&nbsp;&nbsp; Great<br>'
        "<br><em>Dot size also grows with magnitude</em>"
    )
    features.append({
        "type": "Feature",
        "geometry": {"type": "Point", "coordinates": [-165.0, -55.0]},
        "properties": {
            "title":         "Legend",
            "description":   legend_desc,
            "marker-symbol": "information",
            "marker-color":  "#555555",
            "marker-size":   "medium",
        },
    })

    with open(filename, "w") as f:
        json.dump({"type": "FeatureCollection", "features": features}, f)
    print(f"    saved → {filename}  ({len(features)-1:,} earthquakes + 1 legend)")


# =============================================================================
# SECTION 6 – DATE-STAMPED SNAPSHOTS (TODO 5)
# =============================================================================

SNAPSHOT_DIR = "snapshots"   # dated copies land here; *.png/*.geojson stay gitignored


def archive_snapshot(filename: str, date_stamp: str) -> None:
    """
    Copy a freshly-written output file into SNAPSHOT_DIR with the run's UTC
    date stamped onto its name (e.g. map_epicenters_2026-08-15.png).

    The canonical, undated filenames (used by index.html and the README) are
    left untouched so existing links keep working; this just builds up a
    dated history alongside them for tracking change over time.
    """
    if not os.path.exists(filename):
        return
    os.makedirs(SNAPSHOT_DIR, exist_ok=True)
    stem, ext = os.path.splitext(os.path.basename(filename))
    dated_name = f"{stem}_{date_stamp}{ext}"
    shutil.copy2(filename, os.path.join(SNAPSHOT_DIR, dated_name))


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    date_stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    # --- fetch ---
    geojson = fetch_earthquakes(
        MIN_MAGNITUDE, DAYS_BACK,
        min_lat=MIN_LAT, max_lat=MAX_LAT,
        min_lon=MIN_LON, max_lon=MAX_LON,
    )

    # --- parse ---
    df = parse_to_dataframe(geojson)

    # --- largest event ---
    print_largest_event(df)

    # --- plots ---
    # Build a filename prefix so region runs don't overwrite the global PNGs.
    prefix = f"{REGION_NAME.lower().replace(' ', '_')}_" if REGION_NAME else ""
    print("[4/4] Generating plots …", flush=True)

    map_file      = f"{prefix}map_epicenters.png"
    gr_file       = f"{prefix}gutenberg_richter.png"
    gr_stats_file = f"{prefix}gr_stats.json"
    depth_file    = f"{prefix}depth_histogram.png"
    geojson_file  = f"{prefix}earthquakes.geojson"

    plot_world_map(df, filename=map_file)

    # GR fit uses only consistently-typed magnitudes (see MAG_TYPES); the map,
    # depth histogram, aftershock series, and Benioff plot use every event.
    gr_df = filter_by_magnitude_type(df, MAG_TYPES)
    gr_stats = plot_gutenberg_richter(gr_df, filename=gr_file)
    with open(gr_stats_file, "w") as f:
        json.dump(gr_stats, f, indent=2)

    plot_depth_histogram(df, filename=depth_file)

    # --- interactive GeoJSON ---
    save_geojson(df, filename=geojson_file)

    for f in (map_file, gr_file, gr_stats_file, depth_file, geojson_file):
        archive_snapshot(f, date_stamp)

    # --- aftershock time series (TODO 3) ---
    aftershock_file = f"{prefix}aftershock_series.png"
    plot_aftershock_series(df, filename=aftershock_file)
    archive_snapshot(aftershock_file, date_stamp)

    # --- Benioff zone depth gradient (TODO 4) ---
    if BENIOFF_TRANSECT is not None:
        benioff_file = f"{prefix}benioff_zone.png"
        plot_benioff_zone(df, BENIOFF_TRANSECT, filename=benioff_file)
        archive_snapshot(benioff_file, date_stamp)

    # --- b-value comparison (TODO 2) ---
    if COMPARE_REGIONS is not None:
        r_a, r_b = COMPARE_REGIONS
        print("[5/5] Comparing b-values …", flush=True)

        geo_a = fetch_earthquakes(
            MIN_MAGNITUDE, DAYS_BACK,
            min_lat=r_a["min_lat"], max_lat=r_a["max_lat"],
            min_lon=r_a["min_lon"], max_lon=r_a["max_lon"],
        )
        df_a = parse_to_dataframe(geo_a)

        geo_b = fetch_earthquakes(
            MIN_MAGNITUDE, DAYS_BACK,
            min_lat=r_b["min_lat"], max_lat=r_b["max_lat"],
            min_lon=r_b["min_lon"], max_lon=r_b["max_lon"],
        )
        df_b = parse_to_dataframe(geo_b)

        gr_comparison_file = f"{prefix}gr_comparison.png"
        plot_gr_comparison(
            df_a, r_a["name"],
            df_b, r_b["name"],
            filename=gr_comparison_file,
        )
        archive_snapshot(gr_comparison_file, date_stamp)

    print("\nDone!  Push earthquakes.geojson to GitHub to view the interactive map.")


if __name__ == "__main__":
    main()


# =============================================================================
# TODOs – ideas for extending this project
# =============================================================================

# DONE 1 – Geographic bounding box filter
#   Set REGION_NAME / MIN_LAT / MAX_LAT / MIN_LON / MAX_LON in the CONFIG block.

# DONE 2 – Compare b-values between two regions
#   Set COMPARE_REGIONS in the CONFIG block to a list of two region dicts.
#   The script fetches each region separately and overlays both GR fits on
#   gr_comparison.png.  A lower b-value means relatively more large quakes.

# DONE 3 – Aftershock time series
#   Finds the largest event, collects aftershocks within 200 km / 90 days,
#   plots magnitude vs. time, and overlays an Omori–Utsu decay fit on the
#   daily rate — saved to aftershock_series.png.

# DONE 4 – Depth gradient across a subduction zone (Benioff zone)
#   Set BENIOFF_TRANSECT in the CONFIG block to a trench-to-backarc line
#   (start/end lat-lon + corridor half-width).  Projects events onto the
#   line and plots along-track distance vs. depth, fitting a slab dip angle
#   — saved to benioff_zone.png.  Set BENIOFF_TRANSECT = None to skip it.

# DONE 5 – Date-stamped output filenames for tracking change over time
#   Every run copies its outputs into snapshots/ with the UTC run date
#   stamped onto the filename (e.g. snapshots/gutenberg_richter_2026-08-15.png),
#   while the canonical undated files stay in place for index.html/README.

# TODO 6 – Pair catalog events with InSAR coseismic interferograms
#   For M6.5+ events, pull Sentinel-1 SAR pairs (pre/post-event) from the
#   Alaska Satellite Facility (ASF) archive, generate a coseismic
#   interferogram, and compare the deformation fringe pattern's location
#   and inferred fault geometry against the catalog epicenter and depth.
#   Harder than the plots above: needs an ASF/Copernicus data search
#   (asf_search or the ASF Vertex API), an InSAR processing step (e.g.
#   ISCE2, MintPy, or a hosted service — raw Sentinel-1 SLCs are large and
#   the coregistration/unwrapping step is nontrivial), and a way to
#   reconcile the interferogram's phase-derived location/extent with the
#   point-source catalog values.

# TODO 7 – Fit ETAS to the aftershock sequences instead of just plotting them
#   Replace plot_aftershock_series()'s two-parameter Omori-Utsu curve fit
#   with a real Epidemic-Type Aftershock Sequence (ETAS) model: maximize the
#   point-process log-likelihood over (mu, K, alpha, c, p) given the observed
#   event times and magnitudes, instead of log-log regression on a binned
#   daily rate. A fitted ETAS model gives testable forecasts — e.g. expected
#   event counts in a future time window — that can be scored against what
#   the catalog actually recorded, which the current curve fit can't do.
