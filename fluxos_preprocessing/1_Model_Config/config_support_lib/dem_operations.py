#!/usr/bin/env python3
"""
DEM operations for FLUXOS preprocessing.
Handles GeoTIFF reading, downscaling, and ESRI ASCII (.asc) export.
"""

import sys
import numpy as np


def read_geotiff(filepath):
    """
    Read a GeoTIFF DEM file.

    Parameters
    ----------
    filepath : str
        Path to GeoTIFF file.

    Returns
    -------
    elevation : numpy.ndarray
        2D elevation array (nrows x ncols).
    transform : affine.Affine
        Affine transform (pixel to CRS coordinates).
    crs : rasterio.crs.CRS
        Coordinate reference system.
    nodata : float or None
        NODATA value from the file.
    """
    import rasterio

    with rasterio.open(filepath) as src:
        elevation = src.read(1)  # first band
        transform = src.transform
        crs = src.crs
        nodata = src.nodata

        # Validate CRS is projected (meters), not geographic (degrees)
        if crs is not None and crs.is_geographic:
            print("WARNING: DEM is in geographic coordinates (degrees).")
            print("  For accurate slope computation and mesh sizing, use a projected CRS (e.g., UTM).")
            print("  Consider reprojecting with: gdalwarp -t_srs EPSG:xxxxx input.tif output.tif")

    if nodata is None:
        nodata = -9999.0
        print(f"  No NODATA value found in TIF, defaulting to {nodata}")

    return elevation, transform, crs, nodata


def print_tif_info(filepath):
    """
    Print metadata of a GeoTIFF file.

    Parameters
    ----------
    filepath : str
        Path to GeoTIFF file.
    """
    import rasterio

    with rasterio.open(filepath) as src:
        print(f"\n=== GeoTIFF Info: {filepath} ===")
        print(f"  Dimensions:  {src.width} cols x {src.height} rows")
        print(f"  Bands:       {src.count}")
        print(f"  Pixel size:  {abs(src.transform.a):.4f} x {abs(src.transform.e):.4f}")
        print(f"  CRS:         {src.crs}")
        print(f"  NODATA:      {src.nodata}")

        bounds = src.bounds
        print(f"  Bounding box:")
        print(f"    X: {bounds.left:.2f} to {bounds.right:.2f}")
        print(f"    Y: {bounds.bottom:.2f} to {bounds.top:.2f}")

        data = src.read(1)
        valid = data[data != src.nodata] if src.nodata is not None else data.ravel()
        if len(valid) > 0:
            print(f"  Elevation range: {valid.min():.2f} to {valid.max():.2f}")
            print(f"  Valid cells:     {len(valid)} / {data.size} ({100*len(valid)/data.size:.1f}%)")
        else:
            print("  WARNING: No valid elevation data found")

        print()


def downscale_dem(elevation, src_transform, target_resolution, nodata_value):
    """
    Resample DEM to a target resolution using bilinear interpolation.

    Parameters
    ----------
    elevation : numpy.ndarray
        2D elevation array.
    src_transform : affine.Affine
        Source affine transform.
    target_resolution : float
        Target cell size in CRS units (typically meters).
    nodata_value : float
        NODATA value.

    Returns
    -------
    new_elevation : numpy.ndarray
        Resampled elevation array.
    new_transform : affine.Affine
        New affine transform for the resampled grid.
    """
    from scipy.ndimage import zoom
    import rasterio.transform

    src_res = abs(src_transform.a)  # pixel width in CRS units
    zoom_factor = src_res / target_resolution

    nrows, ncols = elevation.shape
    print(f"  Source resolution: {src_res:.4f}")
    print(f"  Target resolution: {target_resolution:.4f}")
    print(f"  Zoom factor: {zoom_factor:.4f}")

    # Create masked copy for interpolation
    elev_float = elevation.astype(np.float64).copy()
    mask = np.isclose(elev_float, nodata_value)
    elev_float[mask] = np.nan

    # Resample using scipy zoom (bilinear interpolation)
    new_elevation = zoom(elev_float, zoom_factor, order=1, mode='nearest')

    # Restore NODATA where NaN
    nodata_mask = np.isnan(new_elevation)
    new_elevation[nodata_mask] = nodata_value

    # Compute new transform
    new_nrows, new_ncols = new_elevation.shape
    xmin = src_transform.c
    ymax = src_transform.f
    xmax = xmin + src_transform.a * ncols
    ymin = ymax + src_transform.e * nrows

    new_transform = rasterio.transform.from_bounds(
        xmin, ymin, xmax, ymax,
        new_ncols, new_nrows
    )

    print(f"  New dimensions: {new_ncols} cols x {new_nrows} rows")

    return new_elevation, new_transform


def write_esri_ascii(filepath, elevation, xllcorner, yllcorner, cellsize, nodata_value):
    """
    Write elevation data in ESRI ASCII raster (.asc) format.

    The format matches what FLUXOS read_geo() expects:
    - 6-line header: NCOLS, NROWS, XLLCORNER, YLLCORNER, CELLSIZE, NODATA_VALUE
    - Data rows from top (north) to bottom (south)

    Parameters
    ----------
    filepath : str
        Output file path.
    elevation : numpy.ndarray
        2D elevation array (nrows x ncols). First row = northernmost.
    xllcorner : float
        X coordinate of lower-left corner.
    yllcorner : float
        Y coordinate of lower-left corner.
    cellsize : float
        Cell size in CRS units.
    nodata_value : float
        NODATA value.
    """
    nrows, ncols = elevation.shape

    with open(filepath, 'w') as f:
        # Header (FLUXOS reads NCOLS/NROWS first, then XLLCORNER, YLLCORNER, CELLSIZE, NODATA)
        f.write(f"NCOLS {ncols}\n")
        f.write(f"NROWS {nrows}\n")
        f.write(f"XLLCORNER {xllcorner:.6f}\n")
        f.write(f"YLLCORNER {yllcorner:.6f}\n")
        f.write(f"CELLSIZE {cellsize:.6f}\n")
        f.write(f"NODATA_VALUE {int(nodata_value)}\n")

        # Data rows: top-to-bottom (north to south), matching standard ESRI ASCII
        for row in range(nrows):
            values = elevation[row, :]
            line = ' '.join(f"{v:.4f}" for v in values)
            f.write(line + '\n')

    print(f"  Written ESRI ASCII: {filepath}")
    print(f"    {ncols} cols x {nrows} rows, cellsize={cellsize:.4f}")


def write_dummy_asc(filepath, xllcorner, yllcorner, cellsize, nodata_value=-9999):
    """
    Write a minimal (10x10) ESRI ASCII file with NODATA values.

    Used for the triangular mesh path where the DEM is embedded in .msh vertex z,
    but main.cpp still requires a valid DEM_FILE for get_domain_size().

    Parameters
    ----------
    filepath : str
        Output file path.
    xllcorner : float
        X coordinate of lower-left corner.
    yllcorner : float
        Y coordinate of lower-left corner.
    cellsize : float
        Cell size.
    nodata_value : float
        NODATA value.
    """
    nrows = ncols = 10
    elevation = np.full((nrows, ncols), nodata_value)

    write_esri_ascii(filepath, elevation, xllcorner, yllcorner, cellsize, nodata_value)
    print(f"  (Dummy ASC for triangular mesh: {nrows}x{ncols})")


def cmd_dem(args):
    """
    CLI handler for the 'dem' subcommand.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    """
    # Info mode: just print metadata and exit
    if args.info:
        print_tif_info(args.input)
        return

    # Must have output specified
    if not args.output_asc:
        print("ERROR: --output-asc is required (unless using --info)")
        sys.exit(1)

    # Read GeoTIFF
    print(f"\nReading GeoTIFF: {args.input}")
    elevation, transform, crs, nodata = read_geotiff(args.input)
    print(f"  Loaded: {elevation.shape[1]} cols x {elevation.shape[0]} rows")

    src_res = abs(transform.a)

    # Downscale if target resolution specified
    if args.resolution:
        if args.resolution <= 0:
            print("ERROR: --resolution must be positive")
            sys.exit(1)

        print(f"\nDownscaling to {args.resolution}m resolution...")
        elevation, transform = downscale_dem(
            elevation, transform, args.resolution, nodata
        )
        cellsize = args.resolution
    else:
        cellsize = src_res
        print(f"  Using native resolution: {cellsize:.4f}")

    # Compute lower-left corner from transform
    nrows, ncols = elevation.shape
    xllcorner = transform.c  # upper-left x
    yllcorner = transform.f + transform.e * nrows  # upper-left y + nrows * (negative pixel height)

    # Write ASC
    print(f"\nWriting ESRI ASCII: {args.output_asc}")
    write_esri_ascii(args.output_asc, elevation, xllcorner, yllcorner, cellsize, nodata)

    print("\nDone!")
