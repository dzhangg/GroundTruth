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
from datetime import datetime, timedelta, timezone

import matplotlib; matplotlib.use("Agg")  # must precede pyplot import
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
            "place":     props["place"],
            "longitude": coords[0],
            "latitude":  coords[1],
            "depth_km":  coords[2],
        })

    df = pd.DataFrame(rows)

    # Convert types so later math works cleanly
    df["magnitude"] = pd.to_numeric(df["magnitude"], errors="coerce")
    df["depth_km"]  = pd.to_numeric(df["depth_km"],  errors="coerce")

    # Drop rows where essential values are missing
    df = df.dropna(subset=["magnitude", "latitude", "longitude", "depth_km"])

    print(f"    → DataFrame has {len(df)} rows × {len(df.columns)} columns")
    return df


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
) -> None:
    """
    Gutenberg-Richter (GR) frequency-magnitude plot.

    The GR relation states:   log10( N(≥M) ) = a − b·M
    where N(≥M) is the number of earthquakes with magnitude ≥ M.
    The 'b-value' (≈ 1 for most regions) describes how many small
    quakes accompany each large one.  We estimate it by fitting a
    straight line to the log10(cumulative count) vs. M curve.
    """
    # Build a range of magnitude bins from the data minimum to maximum
    m_min = math.floor(df["magnitude"].min() * 10) / 10   # round down to 0.1
    m_max = math.ceil(df["magnitude"].max()  * 10) / 10   # round up   to 0.1
    magnitudes = np.arange(m_min, m_max + 0.1, 0.1)       # step = 0.1

    # For each magnitude threshold M, count how many quakes have mag ≥ M
    counts = np.array(
        [(df["magnitude"] >= m).sum() for m in magnitudes], dtype=float
    )

    # Keep only bins with at least one earthquake so log10 doesn't blow up
    valid = counts > 0
    mag_valid   = magnitudes[valid]
    log10_count = np.log10(counts[valid])

    # Fit a degree-1 polynomial (straight line) to the log10 counts vs magnitude.
    # np.polyfit returns [slope, intercept]; slope is −b in the GR relation.
    coeffs = np.polyfit(mag_valid, log10_count, 1)
    b_value = -coeffs[0]          # b-value is the negative of the slope
    a_value =  coeffs[1]

    # Evaluate the fitted line at a smooth set of magnitude values
    fit_x = np.linspace(mag_valid.min(), mag_valid.max(), 200)
    fit_y = np.polyval(coeffs, fit_x)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(mag_valid, log10_count, color="steelblue", s=30,
               label="Observed cumulative count", zorder=3)

    ax.plot(fit_x, fit_y, color="tomato", linewidth=2,
            label=f"GR fit  (b = {b_value:.2f})")

    ax.set_xlabel("Magnitude (M)", fontsize=12)
    ax.set_ylabel("log₁₀ N(≥M)", fontsize=12)
    region_label = f" — {REGION_NAME}" if REGION_NAME else ""
    ax.set_title(
        f"Gutenberg-Richter Plot{region_label} — Last {DAYS_BACK} Days (M ≥ {MIN_MAGNITUDE})",
        fontsize=13,
    )
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"    saved → {filename}")


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


# =============================================================================
# SECTION 5 – INTERACTIVE GEOJSON (GitHub / Azure Maps viewer)
# =============================================================================

def save_geojson(df: pd.DataFrame, filename: str = "earthquakes.geojson") -> None:
    """
    Write the DataFrame to a GeoJSON FeatureCollection file.

    GitHub automatically renders .geojson files as an interactive map using
    Azure Maps (OSM tiles + Leaflet).  Each earthquake becomes a clickable
    point whose popup shows magnitude, place, depth, and date.

    GeoJSON Point coordinates must be [longitude, latitude] — depth is carried
    in 'properties' rather than as a third coordinate so GitHub renders it
    correctly.

    We also add a 'marker-color' property using the standard GeoJSON Simple
    Style spec that GitHub's renderer understands:
      shallow  < 70 km  → yellow  #f5c518
      intermediate 70–300 km → orange #f07800
      deep    > 300 km  → purple  #8b00d4
    """
    def depth_color(depth_km: float) -> str:
        if depth_km < 70:
            return "#f5c518"
        if depth_km < 300:
            return "#f07800"
        return "#8b00d4"

    def mag_symbol_size(mag: float) -> str:
        # GitHub Simple Style supports "small", "medium", "large"
        if mag < 5.5:
            return "small"
        if mag < 7.0:
            return "medium"
        return "large"

    features = []
    for _, row in df.iterrows():
        features.append({
            "type": "Feature",
            "geometry": {
                "type":        "Point",
                "coordinates": [row["longitude"], row["latitude"]],
            },
            "properties": {
                # Data fields shown in the popup
                "magnitude": row["magnitude"],
                "place":     row["place"],
                "depth_km":  round(row["depth_km"], 1),
                "date":      row["time"].strftime("%Y-%m-%d"),
                # GitHub Simple Style fields — control dot colour and size
                "marker-color": depth_color(row["depth_km"]),
                "marker-size":  mag_symbol_size(row["magnitude"]),
            },
        })

    with open(filename, "w") as f:
        json.dump({"type": "FeatureCollection", "features": features}, f)
    print(f"    saved → {filename}  ({len(features):,} features)")


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
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
    plot_world_map(df,            filename=f"{prefix}map_epicenters.png")
    plot_gutenberg_richter(df,    filename=f"{prefix}gutenberg_richter.png")
    plot_depth_histogram(df,      filename=f"{prefix}depth_histogram.png")

    # --- interactive GeoJSON ---
    save_geojson(df, filename=f"{prefix}earthquakes.geojson")

    print("\nDone!  Push earthquakes.geojson to GitHub to view the interactive map.")


if __name__ == "__main__":
    main()


# =============================================================================
# TODOs – ideas for extending this project
# =============================================================================

# DONE 1 – Geographic bounding box filter
#   Set REGION_NAME / MIN_LAT / MAX_LAT / MIN_LON / MAX_LON in the CONFIG block.

# TODO 2 – Compare b-values between two regions
#   Define two bounding boxes (e.g. Japan vs. the Andes), fetch each separately,
#   fit a GR line to each, and overlay both on the same plot.  A lower b-value
#   suggests a region with more large quakes relative to small ones.

# TODO 3 – Aftershock time series
#   After a big mainshock, filter the DataFrame to quakes within ~200 km of the
#   epicentre and within 90 days after.  Plot magnitude vs. time as a scatter
#   plot.  Overlay the Omori–Utsu decay law  n(t) = K / (c + t)^p  fitted to
#   the aftershock rate — this describes how quickly aftershocks die off.
