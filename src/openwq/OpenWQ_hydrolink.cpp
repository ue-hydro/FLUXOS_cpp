// Copyright 2020, Diogo Costa, diogo.pinhodacosta@canada.ca
// This file is part of OpenWQ model.

// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Stub implementations for compiling FLUXOS without the full OpenWQ library.
// When OpenWQ is fully linked, replace this file with the full implementation.

#include "OpenWQ_hydrolink.h"

void openwq_hydrolink::openwq_decl(
    OpenWQ_couplercalls& /*OpenWQ_couplercalls*/,
    OpenWQ_hostModelconfig& /*OpenWQ_hostModelconfig*/,
    OpenWQ_json& /*OpenWQ_json*/,
    OpenWQ_wqconfig& /*OpenWQ_wqconfig*/,
    OpenWQ_units& /*OpenWQ_units*/,
    OpenWQ_utils& /*OpenWQ_utils*/,
    OpenWQ_readjson& /*OpenWQ_readjson*/,
    OpenWQ_vars& /*OpenWQ_vars*/,
    OpenWQ_initiate& /*OpenWQ_initiate*/,
    OpenWQ_watertransp& /*OpenWQ_watertransp*/,
    OpenWQ_chem& /*OpenWQ_chem*/,
    OpenWQ_extwatflux_ss& /*OpenWQ_extwatflux_ss*/,
    OpenWQ_output& /*OpenWQ_output*/,
    std::string /*openwq_masterfile*/,
    unsigned long /*MROWS*/,
    unsigned long /*MCOLS*/)
{
    // Stub: no-op without full OpenWQ library
}

void openwq_hydrolink::openwq_time_start(
    OpenWQ_couplercalls& /*OpenWQ_couplercalls*/,
    OpenWQ_hostModelconfig& /*OpenWQ_hostModelconfig*/,
    OpenWQ_json& /*OpenWQ_json*/,
    OpenWQ_wqconfig& /*OpenWQ_wqconfig*/,
    OpenWQ_units& /*OpenWQ_units*/,
    OpenWQ_utils& /*OpenWQ_utils*/,
    OpenWQ_readjson& /*OpenWQ_readjson*/,
    OpenWQ_vars& /*OpenWQ_vars*/,
    OpenWQ_initiate& /*OpenWQ_initiate*/,
    OpenWQ_watertransp& /*OpenWQ_watertransp*/,
    OpenWQ_chem& /*OpenWQ_chem*/,
    OpenWQ_extwatflux_ss& /*OpenWQ_extwatflux_ss*/,
    OpenWQ_solver& /*OpenWQ_solver*/,
    OpenWQ_output& /*OpenWQ_output*/,
    GlobVar& /*GlobVar_fluxos*/)
{
    // Stub: no-op without full OpenWQ library
}

void openwq_hydrolink::run_space_in(
    GlobVar& /*GlobVar*/,
    OpenWQ_couplercalls& /*OpenWQ_couplercalls*/,
    OpenWQ_hostModelconfig& /*OpenWQ_hostModelconfig*/,
    OpenWQ_json& /*OpenWQ_json*/,
    OpenWQ_wqconfig& /*OpenWQ_wqconfig*/,
    OpenWQ_units& /*OpenWQ_units*/,
    OpenWQ_utils& /*OpenWQ_utils*/,
    OpenWQ_readjson& /*OpenWQ_readjson*/,
    OpenWQ_vars& /*OpenWQ_vars*/,
    OpenWQ_initiate& /*OpenWQ_initiate*/,
    OpenWQ_watertransp& /*OpenWQ_watertransp*/,
    OpenWQ_chem& /*OpenWQ_chem*/,
    OpenWQ_extwatflux_ss& /*OpenWQ_extwatflux_ss*/,
    OpenWQ_solver& /*OpenWQ_solver*/,
    OpenWQ_output& /*OpenWQ_output*/,
    std::string /*source_EWF_name*/,
    int /*ix_r*/, int /*iy_r*/,
    double /*wflux_s2r*/)
{
    // Stub: no-op without full OpenWQ library
}

void openwq_hydrolink::openwq_time_end(
    OpenWQ_couplercalls& /*OpenWQ_couplercalls*/,
    OpenWQ_hostModelconfig& /*OpenWQ_hostModelconfig*/,
    OpenWQ_json& /*OpenWQ_json*/,
    OpenWQ_wqconfig& /*OpenWQ_wqconfig*/,
    OpenWQ_units& /*OpenWQ_units*/,
    OpenWQ_utils& /*OpenWQ_utils*/,
    OpenWQ_readjson& /*OpenWQ_readjson*/,
    OpenWQ_vars& /*OpenWQ_vars*/,
    OpenWQ_initiate& /*OpenWQ_initiate*/,
    OpenWQ_watertransp& /*OpenWQ_watertransp*/,
    OpenWQ_chem& /*OpenWQ_chem*/,
    OpenWQ_extwatflux_ss& /*OpenWQ_extwatflux_ss*/,
    OpenWQ_solver& /*OpenWQ_solver*/,
    OpenWQ_output& /*OpenWQ_output*/,
    GlobVar& /*GlobVar_fluxos*/)
{
    // Stub: no-op without full OpenWQ library
}

time_t openwq_hydrolink::getSimTime(
    OpenWQ_wqconfig& /*OpenWQ_wqconfig*/,
    OpenWQ_units& /*OpenWQ_units*/,
    std::string /*fluxos_sim_start_time_str*/,
    double fluxos_time_secs)
{
    // Stub: return seconds as time_t
    return static_cast<time_t>(fluxos_time_secs);
}
