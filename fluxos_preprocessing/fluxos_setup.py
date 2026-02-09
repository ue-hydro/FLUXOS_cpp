#!/usr/bin/env python3
"""
FLUXOS Preprocessing Tool

Converts GeoTIFF DEMs to FLUXOS-compatible formats and generates configuration files.

Usage:
    python fluxos_setup.py dem --input dem.tif --info
    python fluxos_setup.py dem --input dem.tif --output-asc output.asc --resolution 10
    python fluxos_setup.py mesh --input dem.tif --output domain.msh --min-size 5 --max-size 50
    python fluxos_setup.py config --dem-file output.asc --mesh-type regular --output modset.json
"""

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="fluxos_setup",
        description="FLUXOS preprocessing: DEM conversion, adaptive mesh generation, config creation"
    )
    subparsers = parser.add_subparsers(dest="command", required=True,
                                        help="Available subcommands")

    # =========================================================================
    # dem subcommand
    # =========================================================================
    dem_parser = subparsers.add_parser(
        "dem",
        help="DEM manipulation from GeoTIFF: inspect, downscale, export to ASC"
    )
    dem_parser.add_argument("--input", required=True,
                            help="Input GeoTIFF DEM file")
    dem_parser.add_argument("--output-asc",
                            help="Output ESRI ASCII (.asc) file")
    dem_parser.add_argument("--resolution", type=float,
                            help="Target cell size in meters (for downscaling)")
    dem_parser.add_argument("--info", action="store_true",
                            help="Print GeoTIFF metadata and exit")

    # =========================================================================
    # mesh subcommand
    # =========================================================================
    mesh_parser = subparsers.add_parser(
        "mesh",
        help="Generate adaptive triangular mesh from DEM (slope-based refinement)"
    )
    mesh_parser.add_argument("--input", required=True,
                             help="Input GeoTIFF DEM file")
    mesh_parser.add_argument("--output", default="domain.msh",
                             help="Output Gmsh .msh file (default: domain.msh)")
    mesh_parser.add_argument("--min-size", type=float, default=5.0,
                             help="Minimum triangle edge length in meters (default: 5.0)")
    mesh_parser.add_argument("--max-size", type=float, default=50.0,
                             help="Maximum triangle edge length in meters (default: 50.0)")
    mesh_parser.add_argument("--slope-factor", type=float, default=1.0,
                             help="Slope-to-refinement aggressiveness: "
                                  "0 = uniform mesh, 1 = full adaptation (default: 1.0)")
    mesh_parser.add_argument("--boundary-tags", default="wall",
                             help="Comma-separated boundary tags (default: wall)")

    # =========================================================================
    # config subcommand
    # =========================================================================
    config_parser = subparsers.add_parser(
        "config",
        help="Generate FLUXOS JSON configuration file (modset.json)"
    )
    config_parser.add_argument("--dem-file", required=True,
                               help="Path to DEM file (.asc) for config")
    config_parser.add_argument("--mesh-type",
                               choices=["regular", "triangular"], default="regular",
                               help="Mesh type (default: regular)")
    config_parser.add_argument("--mesh-file",
                               help="Path to mesh file (.msh) - required for triangular")
    config_parser.add_argument("--mesh-format", default="gmsh",
                               help="Mesh file format (default: gmsh)")
    config_parser.add_argument("--meteo-file",
                               help="Path to meteo forcing file (.fluxos)")
    config_parser.add_argument("--inflow-file",
                               help="Path to inflow forcing file (.fluxos)")
    config_parser.add_argument("--inflow-x", type=float,
                               help="Inflow discharge X coordinate")
    config_parser.add_argument("--inflow-y", type=float,
                               help="Inflow discharge Y coordinate")
    config_parser.add_argument("--sim-start", default="2000-01-01 00:00:00",
                               help="Simulation start datetime (default: 2000-01-01 00:00:00)")
    config_parser.add_argument("--output-step", type=int, default=3600,
                               help="Output print step in seconds (default: 3600)")
    config_parser.add_argument("--output-folder", default="Results/",
                               help="Output folder for results (default: Results/)")
    config_parser.add_argument("--roughness", type=float, default=0.005,
                               help="Roughness height in meters (default: 0.005)")
    config_parser.add_argument("--h-min-print", type=float, default=0.005,
                               help="Minimum water depth to print (default: 0.005)")
    config_parser.add_argument("--restart", action="store_true",
                               help="Enable restart from previous simulation")
    config_parser.add_argument("--ade-status", action="store_true",
                               help="Enable ADE transport module")
    config_parser.add_argument("--ade-dcoef", type=float, default=0.5,
                               help="ADE dispersion coefficient (default: 0.5)")
    config_parser.add_argument("--openwq-status", action="store_true",
                               help="Enable OpenWQ module")
    config_parser.add_argument("--openwq-masterfile", default="",
                               help="OpenWQ master config file path")
    config_parser.add_argument("--boundary-conditions",
                               help="Boundary conditions as JSON string")
    config_parser.add_argument("--output", default="modset.json",
                               help="Output JSON config file (default: modset.json)")

    # =========================================================================
    # Parse and dispatch
    # =========================================================================
    args = parser.parse_args()

    if args.command == "dem":
        from dem_operations import cmd_dem
        cmd_dem(args)

    elif args.command == "mesh":
        from mesh_generation import cmd_mesh
        cmd_mesh(args)

    elif args.command == "config":
        from config_generator import cmd_config
        cmd_config(args)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
