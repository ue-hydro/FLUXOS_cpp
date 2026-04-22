# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of the FLUXOS model.

# This program, FLUXOS, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
On-demand DEM fetcher for the FLUXOS model-config template.

Two backends are supported:

1. **OpenTopography global-DEM API** — single HTTPS endpoint serving ~10
   global products (SRTM 30/90 m, Copernicus GLO-30/90, ALOS AW3D-30,
   NASADEM, EU-DTM, GEBCO, …). Requires a free API key from
   <https://portal.opentopography.org/>.

2. **USGS 3DEP** (US only) — via the ``py3dep`` package (optional dependency).
   Provides 10 m nationwide, 3 m or 1 m where LiDAR coverage exists — the
   highest-resolution option for US flood modelling.

Results are cached locally keyed by ``(provider, bbox)`` so repeated runs of
the template do not re-hit the network.

The fetched GeoTIFF is reprojected to a projected CRS (default auto-UTM
based on the bbox centre) before being handed off to the regular
preprocessing pipeline. FLUXOS's slope / mesh code requires metric units.
"""

from __future__ import annotations

import hashlib
import json
import os
import urllib.error
import urllib.request
from typing import Tuple

# ---------------------------------------------------------------------------
#  Provider registry
# ---------------------------------------------------------------------------

# Each OpenTopography entry maps to the ``demtype`` parameter in their API.
# See https://portal.opentopography.org/apidocs/
OPENTOPO_PROVIDERS: dict[str, dict] = {
    "SRTMGL1": {"label": "SRTM 30 m (NASA)",            "native_res_m": 30,
                "vertical_rmse_m": 16, "coverage": "global (lat ±60°)"},
    "SRTMGL3": {"label": "SRTM 90 m (NASA)",            "native_res_m": 90,
                "vertical_rmse_m": 16, "coverage": "global (lat ±60°)"},
    "COP30":   {"label": "Copernicus GLO-30 (ESA)",     "native_res_m": 30,
                "vertical_rmse_m": 4,  "coverage": "global"},
    "COP90":   {"label": "Copernicus GLO-90 (ESA)",     "native_res_m": 90,
                "vertical_rmse_m": 4,  "coverage": "global"},
    "AW3D30":  {"label": "ALOS AW3D30 (JAXA)",          "native_res_m": 30,
                "vertical_rmse_m": 5,  "coverage": "global"},
    "NASADEM": {"label": "NASADEM (NASA, SRTM-reproc.)", "native_res_m": 30,
                "vertical_rmse_m": 10, "coverage": "global (lat ±60°)"},
    "EU_DTM":  {"label": "EU-DTM 10 m (Copernicus LS)", "native_res_m": 10,
                "vertical_rmse_m": 3,  "coverage": "Europe (EEA39)"},
    "GEBCOIceTopo": {"label": "GEBCO ice-surface",      "native_res_m": 463,
                "vertical_rmse_m": None, "coverage": "global ocean + land"},
}

# USGS 3DEP products (US only) — resolved at download time via py3dep.
USGS_PROVIDERS: dict[str, dict] = {
    "USGS_10M": {"label": "USGS 3DEP 10 m (1/3 arc-sec)", "native_res_m": 10,
                 "vertical_rmse_m": 1, "coverage": "CONUS"},
    "USGS_3M":  {"label": "USGS 3DEP 3 m (1/9 arc-sec)",  "native_res_m": 3,
                 "vertical_rmse_m": 0.5, "coverage": "CONUS (partial)"},
    "USGS_1M":  {"label": "USGS 3DEP 1 m LiDAR",          "native_res_m": 1,
                 "vertical_rmse_m": 0.15, "coverage": "US LiDAR-covered only"},
}

ALL_PROVIDERS: dict[str, dict] = {**OPENTOPO_PROVIDERS, **USGS_PROVIDERS}


def provider_info(name: str) -> dict:
    if name not in ALL_PROVIDERS:
        raise ValueError(f"Unknown DEM provider: {name!r}. "
                         f"Valid: {sorted(ALL_PROVIDERS)}")
    return {**ALL_PROVIDERS[name], "name": name}


# ---------------------------------------------------------------------------
#  Cache key
# ---------------------------------------------------------------------------

def _cache_filename(provider: str, bbox: Tuple[float, float, float, float]) -> str:
    key = f"{provider}|{bbox[0]:.6f}|{bbox[1]:.6f}|{bbox[2]:.6f}|{bbox[3]:.6f}"
    short = hashlib.sha1(key.encode()).hexdigest()[:12]
    return f"{provider}_{short}.tif"


# ---------------------------------------------------------------------------
#  Fetchers
# ---------------------------------------------------------------------------

_OPENTOPO_URL = "https://portal.opentopography.org/API/globaldem"


def _fetch_opentopography(provider: str, bbox_wgs84: Tuple[float, float, float, float],
                          out_path: str, api_key: str) -> str:
    lon_min, lat_min, lon_max, lat_max = bbox_wgs84
    if not api_key:
        raise RuntimeError(
            "OpenTopography API key missing. Get a free key at "
            "https://portal.opentopography.org/ and expose it via the "
            "environment variable named in `download_api_key_env` "
            "(default: OPENTOPOGRAPHY_API_KEY).")

    url = (f"{_OPENTOPO_URL}?demtype={provider}"
           f"&south={lat_min}&north={lat_max}"
           f"&west={lon_min}&east={lon_max}"
           f"&outputFormat=GTiff&API_Key={api_key}")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    print(f"      fetching {provider} for bbox {bbox_wgs84} …")
    try:
        urllib.request.urlretrieve(url, out_path)
    except urllib.error.HTTPError as e:
        # The API returns plain-text error messages in the body, which the
        # default urllib exception hides — surface them to the user.
        body = e.read().decode("utf-8", errors="replace")[:500]
        raise RuntimeError(
            f"OpenTopography HTTP {e.code}: {body}\n"
            f"URL: {url.split('&API_Key=')[0]}&API_Key=***") from None
    if os.path.getsize(out_path) < 1024:
        # Response that isn't a GeoTIFF (e.g., HTML error page)
        with open(out_path) as f:
            body = f.read(500)
        os.remove(out_path)
        raise RuntimeError(f"OpenTopography returned a non-raster response: {body}")
    return out_path


def _fetch_usgs_3dep(provider: str, bbox_wgs84: Tuple[float, float, float, float],
                     out_path: str) -> str:
    """Fetch 10 m / 3 m / 1 m DEM from USGS 3DEP via py3dep."""
    try:
        import py3dep  # type: ignore
    except ImportError as e:
        raise ImportError(
            "py3dep is required for USGS 3DEP downloads. "
            "Install with `pip install py3dep`.") from e

    # py3dep accepts (minx, miny, maxx, maxy) in WGS84
    res_map = {"USGS_10M": 10, "USGS_3M": 3, "USGS_1M": 1}
    res = res_map[provider]
    print(f"      fetching USGS 3DEP {res} m for bbox {bbox_wgs84} …")
    da = py3dep.get_dem(bbox_wgs84, res)   # xarray.DataArray in WGS84
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    da.rio.to_raster(out_path)
    return out_path


# ---------------------------------------------------------------------------
#  UTM reprojection
# ---------------------------------------------------------------------------

def _auto_utm_epsg(bbox_wgs84: Tuple[float, float, float, float]) -> int:
    lon_min, lat_min, lon_max, lat_max = bbox_wgs84
    cx = 0.5 * (lon_min + lon_max)
    cy = 0.5 * (lat_min + lat_max)
    zone = int((cx + 180) // 6) + 1
    return (32600 if cy >= 0 else 32700) + zone


def reproject_to_projected(src_path: str, dst_path: str,
                           target_crs: str | int = "auto") -> str:
    """
    Reproject a GeoTIFF to a projected metric CRS (default auto-UTM).
    Uses bilinear resampling which is appropriate for continuous elevation.
    """
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    with rasterio.open(src_path) as src:
        if target_crs == "auto":
            bounds = src.bounds
            if src.crs is not None and src.crs.is_geographic:
                bbox = (bounds.left, bounds.bottom, bounds.right, bounds.top)
            else:
                # Already projected — just pass through
                if src_path != dst_path:
                    with open(src_path, "rb") as fin, open(dst_path, "wb") as fout:
                        fout.write(fin.read())
                return dst_path
            epsg = _auto_utm_epsg(bbox)
            dst_crs = f"EPSG:{epsg}"
        else:
            dst_crs = target_crs if isinstance(target_crs, str) else f"EPSG:{target_crs}"

        if src.crs is not None and str(src.crs) == dst_crs:
            # Same CRS → no need to warp
            if src_path != dst_path:
                with open(src_path, "rb") as fin, open(dst_path, "wb") as fout:
                    fout.write(fin.read())
            return dst_path

        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds
        )
        meta = src.meta.copy()
        meta.update(crs=dst_crs, transform=transform, width=width, height=height)
        with rasterio.open(dst_path, "w", **meta) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.bilinear,
                )
        return dst_path


# ---------------------------------------------------------------------------
#  Public entry
# ---------------------------------------------------------------------------

def fetch_dem(provider: str,
              bbox_wgs84: Tuple[float, float, float, float],
              cache_dir: str,
              api_key: str | None = None,
              target_crs: str | int = "auto",
              force: bool = False) -> dict:
    """
    Download (or return cached) DEM for the given bbox and reproject to a
    projected CRS. Returns a metadata dict with the final GeoTIFF path and
    provider info for the report.

    Parameters
    ----------
    provider : str
        Key from ``OPENTOPO_PROVIDERS`` or ``USGS_PROVIDERS``.
    bbox_wgs84 : (lon_min, lat_min, lon_max, lat_max)
    cache_dir : local dir where raw + reprojected GeoTIFFs are kept.
    api_key : OpenTopography API key (ignored for USGS).
    target_crs : ``"auto"`` (UTM from bbox centre), an EPSG code, or a PROJ string.
    force : bypass the cache.
    """
    info = provider_info(provider)
    os.makedirs(cache_dir, exist_ok=True)

    raw_name = _cache_filename(provider, bbox_wgs84)
    raw_path = os.path.join(cache_dir, "raw_" + raw_name)
    proj_path = os.path.join(cache_dir, raw_name)

    if not os.path.exists(raw_path) or force:
        if provider in OPENTOPO_PROVIDERS:
            _fetch_opentopography(provider, bbox_wgs84, raw_path, api_key or "")
        elif provider in USGS_PROVIDERS:
            _fetch_usgs_3dep(provider, bbox_wgs84, raw_path)
        else:
            raise ValueError(f"Unknown provider: {provider}")
    else:
        print(f"      using cached: {raw_path}")

    reproject_to_projected(raw_path, proj_path, target_crs=target_crs)

    # Read final metadata for the report
    import rasterio
    with rasterio.open(proj_path) as ds:
        final_crs = str(ds.crs)
        w, h = ds.width, ds.height
        final_res = abs(ds.transform.a)

    return {
        **info,
        "provider": provider,
        "bbox_wgs84": list(bbox_wgs84),
        "raw_path": raw_path,
        "output_path": proj_path,
        "final_crs": final_crs,
        "final_resolution_m": float(final_res),
        "final_grid": [int(w), int(h)],
    }
