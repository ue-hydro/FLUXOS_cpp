

// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef OPENWQ_HYDROLINK_INCLUDED
#define OPENWQ_HYDROLINK_INCLUDED

#include "../openwq/openwq/src/OpenWQ_initiate.h"
#include "../openwq/openwq/src/OpenWQ_chem.h"
#include "../openwq/openwq/src/OpenWQ_watertransp.h"
#include "../openwq/openwq/src/OpenWQ_extwatflux_ss.h"
#include "../openwq/openwq/src/OpenWQ_output.h"
#include "../openwq/openwq/src/OpenWQ_units.h"
#include "../openwq/openwq/src/OpenWQ_utils.h"
#include "../openwq/openwq/src/OpenWQ_solver.h"
#include "../openwq/openwq/src/OpenWQ_couplercalls.h"
#include "../openwq/openwq/src/OpenWQ_global.h"
#include "../openwq/openwq/src/OpenWQ_readjson.h"

#include "../fluxos/GlobVar.h"
//GlobVar GlobVar;

class openwq_hydrolink
{

    public:

    void openwq_decl(
        OpenWQ_couplercalls& OpenWQ_couplercalls,
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,
        OpenWQ_readjson& OpenWQ_readjson,            // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
        OpenWQ_watertransp& OpenWQ_watertransp,      // transport modules
        OpenWQ_chem& OpenWQ_chem,                   // biochemistry modules
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
        OpenWQ_output& OpenWQ_output,                // output modules (needed for console/logfile)
        std::string openwq_masterfile,
        unsigned long NROWS,
        unsigned long NCOLS);

    
    void openwq_time_start(
        OpenWQ_couplercalls& OpenWQ_couplercalls,
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,
        OpenWQ_readjson& OpenWQ_readjson,            // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
        OpenWQ_watertransp& OpenWQ_watertransp,      // transport modules
        OpenWQ_chem& OpenWQ_chem,                    // biochemistry modules
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
        OpenWQ_solver& OpenWQ_solver,                // solver module
        OpenWQ_output& OpenWQ_output,                // output modules
        GlobVar& GlobVar);                                     // fluxos h

    void run_space_in(
        GlobVar& GlobVar,
        OpenWQ_couplercalls& OpenWQ_couplercalls,
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,
        OpenWQ_readjson& OpenWQ_readjson,            // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
        OpenWQ_watertransp& OpenWQ_watertransp,      // transport modules
        OpenWQ_chem& OpenWQ_chem,                    // biochemistry modules
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,  // sink and source modules)
        OpenWQ_solver& OpenWQ_solver,                // solver module
        OpenWQ_output& OpenWQ_output,                // output modules
        std::string source_EWF_name,
        int ix_r, int iy_r,
        double wflux_s2r);

    void openwq_time_end(
        OpenWQ_couplercalls& OpenWQ_couplercalls,     // Class with all call from coupler
        OpenWQ_hostModelconfig& OpenWQ_hostModelconfig,
        OpenWQ_json& OpenWQ_json,                    // create OpenWQ_json object
        OpenWQ_wqconfig& OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
        OpenWQ_units& OpenWQ_units,                  // functions for unit conversion
        OpenWQ_utils& OpenWQ_utils,
        OpenWQ_readjson& OpenWQ_readjson,            // read json files
        OpenWQ_vars& OpenWQ_vars,
        OpenWQ_initiate& OpenWQ_initiate,            // initiate modules
        OpenWQ_watertransp& OpenWQ_watertransp,      // transport modules
        OpenWQ_chem& OpenWQ_chem,                   // biochemistry modules
        OpenWQ_extwatflux_ss& OpenWQ_extwatflux_ss,        // sink and source modules)
        OpenWQ_solver& OpenWQ_solver,
        OpenWQ_output& OpenWQ_output,
        GlobVar& GlobVar_fluxos);

    time_t getSimTime(
        OpenWQ_wqconfig& OpenWQ_wqconfig,
        OpenWQ_units& OpenWQ_units,
        std::string fluxos_sim_start_time_str, 
        double fluxos_time_secs);


};


#endif