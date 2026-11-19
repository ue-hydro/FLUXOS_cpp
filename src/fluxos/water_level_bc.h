// Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
// This file is part of the FLUXOS model.

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// ---------------------------------------------------------------------------
// Time-varying water-level (stage) boundary condition — REGULAR mesh.
//
// Modset:
//   "BOUNDARY_CONDITIONS": {
//     "WATER_LEVEL": {
//       "STATUS":   true,
//       "FILENAME": "sado_mouth_level.csv",   // t_seconds,level_m (header ok)
//       "EDGE":     "west"                    // see below
//       "VELOCITY": "zero_gradient"           // optional, or "zero"
//     }
//   }
//
// EDGE is one of
//   "west" | "east" | "south" | "north"          whole edge
//   ["west", "south"]                            several whole edges
//   [{"EDGE":"west","FROM":10,"TO":40}, ...]     segments, 1-based cell index
//        along the edge (rows counted from the south for west/east edges,
//        columns counted from the west for south/north edges), inclusive
//   [{"EDGE":"west","FROM_COORD":4.27e6,"TO_COORD":4.28e6}]  segments in map
//        coordinates (northing for west/east, easting for south/north)
//
// Grid orientation (see read_geo): irow = 1..NROWS runs south→north,
// icol = 1..NCOLS runs west→east.
//
// Each time step, after the inflow/meteo forcing and before the flow
// solver, the BC cells get h = max(level(t) − zb, 0) (level linearly
// interpolated in time, held constant outside the series). Momentum is
// copied from the first interior cell as a velocity (zero-gradient), or
// set to zero with VELOCITY = "zero". The legacy "weir" edge-flux override
// is switched off on BC cells so they exchange water with the interior via
// the normal Roe/hydrostatic fluxes; the outer face of a BC cell carries no
// flux. The volume the BC adds/removes is accumulated (wl_bc_vol_*) and
// written to <OUTPUT_FOLDER>/water_level_bc_balance.csv at each print step.
// ---------------------------------------------------------------------------

#ifndef WATER_LEVEL_BC_H_INCLUDED
#define WATER_LEVEL_BC_H_INCLUDED

#include <fstream>
#include "GlobVar.h"

// Parse BOUNDARY_CONDITIONS.WATER_LEVEL, load the series and build the
// BC cell list. Must be called after read_geo (needs NODATA mask and
// geo-referencing). Returns true on error.
bool read_water_level_bc(GlobVar& ds, std::ofstream& logFLUXOSfile);

// Linear interpolation of the level series at time t.
double water_level_at(const GlobVar& ds, double t);

// Impose level(ds.tim) on the BC cells (no-op when disabled).
void apply_water_level_bc(GlobVar& ds, int nchem);

// Append one line (t, level, cumulative volumes) to the balance CSV.
void write_water_level_bc_balance(const GlobVar& ds, double t_print);

#endif
