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
import sys
from datetime import datetime, timedelta, timezone

import matplotlib
matplotlib.use("Agg")          # non-interactive backend — saves PNGs without a display
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import requests

# =============================================================================
# CONFIG  ← change these two numbers to adjust the query
# =============================================================================
MIN_MAGNITUDE = 4.5    # only fetch quakes at or above this magnitude
DAYS_BACK     = 365    # how many days of history to fetch


# =============================================================================
# SECTION 1 – BUILD THE API URL AND DOWNLOAD DATA
# =============================================================================

def fetch_earthquakes(min_mag: float, days_back: int) -> dict:
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

    base_url = "https://earthquake.usgs.gov/fdsnws/event/1/query"

    print(f"[1/4] Fetching M{min_mag}+ earthquakes from {start_time.date()} to "
          f"{end_time.date()} …", flush=True)

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

def plot_world_map(df: pd.DataFrame, filename: str = "map_epicenters.png") -> None:
    """
    Scatter plot of earthquake epicenters on a lon/lat grid.

    Colour encodes depth (shallow = yellow, deep = purple).
    Point size encodes magnitude (bigger quake → bigger dot).
    """
    fig, ax = plt.subplots(figsize=(14, 7))

    # Draw a simple rectangle to represent the globe outline
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_facecolor("#c8e6f5")   # light blue = ocean

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
    ax.set_title(
        f"M{MIN_MAGNITUDE}+ Earthquake Epicenters — Last {DAYS_BACK} Days\n"
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
    ax.set_title(
        f"Gutenberg-Richter Plot — Last {DAYS_BACK} Days (M ≥ {MIN_MAGNITUDE})",
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
    ax.set_title(
        f"Earthquake Depth Distribution — Last {DAYS_BACK} Days (M ≥ {MIN_MAGNITUDE})",
        fontsize=13,
    )
    ax.legend(fontsize=10)
    ax.grid(True, axis="y", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"    saved → {filename}")


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

def main():
    # --- fetch ---
    geojson = fetch_earthquakes(MIN_MAGNITUDE, DAYS_BACK)

    # --- parse ---
    df = parse_to_dataframe(geojson)

    # --- largest event ---
    print_largest_event(df)

    # --- plots ---
    print("[4/4] Generating plots …", flush=True)
    plot_world_map(df)
    plot_gutenberg_richter(df)
    plot_depth_histogram(df)

    print("\nDone!  Three PNG files saved in the current directory.")


if __name__ == "__main__":
    main()


# =============================================================================
# TODOs – ideas for extending this project
# =============================================================================

# TODO 1 – Filter to a geographic bounding box
#   Add config variables MIN_LAT, MAX_LAT, MIN_LON, MAX_LON and pass them to
#   the API as minlatitude / maxlatitude / minlongitude / maxlongitude params.
#   Then re-run the map and GR plot just for that region (e.g. Japan, California).

# TODO 2 – Compare b-values between two regions
#   Define two bounding boxes (e.g. Japan vs. the Andes), fetch each separately,
#   fit a GR line to each, and overlay both on the same plot.  A lower b-value
#   suggests a region with more large quakes relative to small ones.

# TODO 3 – Aftershock time series
#   After a big mainshock, filter the DataFrame to quakes within ~200 km of the
#   epicentre and within 90 days after.  Plot magnitude vs. time as a scatter
#   plot.  Overlay the Omori–Utsu decay law  n(t) = K / (c + t)^p  fitted to
#   the aftershock rate — this describes how quickly aftershocks die off.
