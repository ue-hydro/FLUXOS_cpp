
// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019, Diogo Costa

// This program, FLUXOS, is free software: you can redistribute it and/or modify
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

#include <vector> 
#include <dirent.h> 
#include <memory> 
#include <string.h>
#include <armadillo>
#include <iostream>

#include "GlobVar.h"
#include "common.h"

std::string SplitFilename (
    const std::string& str)
{
  size_t found;
  found=str.find_last_of("/\\");
  return str.substr(0,found);

}

int getIntNumberFromString(
    std::string s){
   std::stringstream str_strm;
   str_strm << s; //convert the string s into stringstream
   std::string temp_str;
   int temp_int;
   while(!str_strm.eof()) {
      str_strm >> temp_str; //take words into temp_str one by one
      if(std::stringstream(temp_str) >> temp_int) { //try to convert string to int
         return temp_int;
      }
      temp_str = ""; //clear temp string
   }
}

double getFloatNumberFromString(
    std::string s){
   std::stringstream str_strm;
   str_strm << s; //convert the string s into stringstream
   std::string temp_str;
   double temp_double;
   while(!str_strm.eof()) {
      str_strm >> temp_str; //take words into temp_str one by one
      if(std::stringstream(temp_str) >> temp_double) { //try to convert string to int
         return temp_double;
      }
      temp_str = ""; //clear temp string
   }
}

// read file names in Results directory
int findLastStep(
    const char *path) 
{

   struct dirent *entry;
   int i, timestart, filenum = 0, simnum;
   std::vector<char*> filenames; //stringvec filenames, filename_i;
   const char *filename_i;
   char *simnum_str_i;
   DIR *dir = opendir(path);
   
   if (dir != NULL) {
        while ((entry = readdir(dir)) != NULL) {
        filenames.push_back(entry->d_name); // storing the file names
        filenum = filenum + 1;
        }
   }
   closedir(dir);
   
   timestart = 0;
   for(i=2;i<filenum;i++){
       filename_i = filenames[i]; //.assign(filenames[i]); //strcpy(filename_i,(char *)(&filenames[i]));
        simnum_str_i = (char*) malloc(sizeof(filename_i)-2);
        strncpy (simnum_str_i, filename_i, sizeof(filename_i)-2);
        simnum = atoi(simnum_str_i);
        timestart = std::fmax(timestart,simnum);
        free(simnum_str_i);
   }
   
   if (filenum >0)
        free(entry); 
   //free(filename_i);
   
   //free(dir);
   //std::cout << "Start time (s): " << timestart << " (initial conditions available)" << std::endl;
   return timestart;
   
}


// get size of the domain
bool get_domain_size(unsigned int *rown, 
    unsigned int *coln, 
    json master_MODSET_local,
    std::ofstream& logFLUXOSfile)
{

    // Local variables
    std::string str, dem_file_temp, msg;
    bool errflag = false;
    
    // Get DEM filepath
    dem_file_temp = master_MODSET_local["DEM_FILE"];

    // Opening the DEM file.
    unsigned int linei,numi;
    std::string line, stri; 
    std::ifstream myfile (dem_file_temp); //opening the file.

    std::string str_nrows = "NROWS";
    std::string str_ncols = "NCOLS";

    if (myfile.is_open()) //if the file is open
    {

        for(linei=1;linei<=2;linei++) // Just need to read the first 2 lines
        {
            getline (myfile,line); //get one line from the file
            stri = line.substr(0,5); 

            // converts string to uppercase for comparison with str_nrows and str_ncols
            std::for_each(stri.begin(), stri.end(), [](char & c) {c = ::toupper(c);});

            numi = getIntNumberFromString(line);

            if (stri.compare(str_nrows) == 0){
                *rown = numi;
            }else if(stri.compare(str_ncols) == 0){
                *coln = numi;
            }
        }

        myfile.close(); //closing the file
        
    }
    else {
        std::cout << "Unable to open file: " + dem_file_temp << std::endl; //if the file is not open output
        errflag = true;
    }

    return errflag;
    
}

// Add meteo at instant t
bool add_meteo(
    GlobVar& ds,
    openwq_hydrolink& openwq_hydrolink,
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
    OpenWQ_output& OpenWQ_output,
    int nchem){

    unsigned int a,irow, icol,meteo_rowi;
    double meteoi,hp;
    double meteo_conci[nchem];
    bool errflag = false;
    std::string openwq_source_EWF_name = "METEO";

    // Return if no METEO_FILE provided
    if (ds.meteo_file.empty())
        return errflag;


    try{
        for (a=0;a<=(*ds.meteo).col(0).n_elem;a++){
            meteo_rowi = a;
            if ((*ds.meteo).at(a,0) > ds.tim){       
                break;
            }
        }
        
        // Convert mm/day to m/s
        meteoi = (*ds.meteo).at(meteo_rowi,1)/(1000.*3600.*24.)*ds.dtfl;

        // Get chem data
        for (int ichem=0;ichem<nchem;ichem++){
            meteo_conci[ichem] = (*ds.meteo).at(meteo_rowi,ichem+2);
        }
        
        for(icol=1;icol<=ds.NCOLS;icol++)
        {
            for(irow=1;irow<=ds.NROWS;irow++)
            {
                if ((*ds.zb).at(irow,icol) != ds.NODATA_VALUE)
                {
                    hp = std::fmax((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0f); // adesolver hp before adding snowmelt  
                    (*ds.z).at(irow,icol) = (*ds.z).at(irow,icol) + meteoi;   
                    (*ds.h)(irow,icol)=std::fmax((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0f);
                    if ((*ds.h)(irow,icol) <= ds.hdry)
                    {
                        (*ds.ldry).at(irow,icol)=0.0f;
                    }
                    (*ds.h0)(irow,icol) = (*ds.h)(irow,icol);

                    // Calc mass balance for all chemcicals
                    if (ds.ade_solver == true && hp!=0.0f)
                    {   
                        for (int ichem=0;ichem<nchem;ichem++){       
                        
                            if (ds.openwq == false){
                                (*ds.conc_SW)[ichem](irow,icol)=((*ds.conc_SW)[ichem](irow,icol)*hp
                                                        + (meteoi * meteo_conci[ichem])
                                                        )/((*ds.h)(irow,icol)); //adesolver (adjustment for snowmelt)   
                            }else{

                                // call openwq adv_in
                                openwq_hydrolink.run_space_in(
                                    ds,
                                    OpenWQ_couplercalls,
                                    OpenWQ_hostModelconfig,
                                    OpenWQ_json,                    // create OpenWQ_json object
                                    OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
                                    OpenWQ_units,                  // functions for unit conversion
                                    OpenWQ_utils,
                                    OpenWQ_readjson,            // read json files
                                    OpenWQ_vars,
                                    OpenWQ_initiate,            // initiate modules
                                    OpenWQ_watertransp,      // transport modules
                                    OpenWQ_chem,                    // biochemistry modules
                                    OpenWQ_extwatflux_ss,  // sink and source modules)
                                    OpenWQ_solver,                // solver module
                                    OpenWQ_output,
                                    openwq_source_EWF_name,
                                    irow, icol,
                                    meteoi);

                            }
                        
                        }    
                    }
                }
            }
        }
        
    }catch (int e){

        std::cout << "problem in 'add_meteo' module" << std::endl; // err message
        errflag = true;
    }

    return errflag;

}


// Add inflow at instant t
bool add_inflow(
    GlobVar& ds,
    openwq_hydrolink& openwq_hydrolink,
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
    OpenWQ_output& OpenWQ_output,
    int nchem){

unsigned int a,irow, icol,inflow_rowi;
double inflow_conci[nchem];
double inflowi,hp;
bool errflag = false;
std::string openwq_source_EWF_name = "INFLOW"; // for openwq ewf

// Return if no METEO_FILE provided
if (ds.inflow_file.empty())
    return errflag;

try{
    for (a=0;a<=(*ds.inflow).col(0).n_elem;a++){
        inflow_rowi = a;
        if ((*ds.inflow).at(a,0) > ds.tim){       
            break;
        }
    }

    irow = ds.inflow_nrow;
    icol = ds.inflow_ncol;
    
    inflowi = (*ds.inflow).at(inflow_rowi,1)*ds.dtfl/(std::pow(ds.dxy,2)); // added as m3/s
    
    // Get chem data
    for (int ichem=0;ichem<nchem;ichem++){
        inflow_conci[ichem] = (*ds.inflow).at(inflow_rowi,ichem+2);
    }

    if ((*ds.zb).at(irow,icol) != ds.NODATA_VALUE)
    {
        if ((*ds.z)(irow,icol) <= (*ds.zb).at(irow,icol)){
            (*ds.z).at(irow,icol) = (*ds.zb).at(irow,icol);
        }
    
        hp = std::fmax((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0f); // adesolver hp before adding snowmelt  
        (*ds.z).at(irow,icol) = (*ds.z).at(irow,icol) + inflowi;   
        (*ds.h)(irow,icol)=std::fmax((*ds.z).at(irow,icol)-(*ds.zb).at(irow,icol),0.0f);

        if ((*ds.h)(irow,icol) <= ds.hdry)
            (*ds.ldry).at(irow,icol) = 1.0f;
        else
            (*ds.ldry).at(irow,icol)= 0.0f;
        
        (*ds.h0)(irow,icol) = (*ds.h)(irow,icol);

        // Calc mass balance for all chemcicals
        if (ds.ade_solver == true && hp!=0.0f)
        {         
            for (int ichem=0;ichem<nchem;ichem++){  

                if (ds.openwq == false){
                    (*ds.conc_SW)[ichem](irow,icol)=((*ds.conc_SW)[ichem](irow,icol)*hp
                                                + (inflowi * inflow_conci[ichem])
                                            )/((*ds.h)(irow,icol)); //adesolver (adjustment for snowmelt)  
                }else{

                    // call openwq adv_in
                    openwq_hydrolink.run_space_in(
                        ds,
                        OpenWQ_couplercalls,
                        OpenWQ_hostModelconfig,
                        OpenWQ_json,                    // create OpenWQ_json object
                        OpenWQ_wqconfig,            // create OpenWQ_wqconfig object
                        OpenWQ_units,                  // functions for unit conversion
                        OpenWQ_utils,
                        OpenWQ_readjson,            // read json files
                        OpenWQ_vars,
                        OpenWQ_initiate,            // initiate modules
                        OpenWQ_watertransp,      // transport modules
                        OpenWQ_chem,                    // biochemistry modules
                        OpenWQ_extwatflux_ss,  // sink and source modules)
                        OpenWQ_solver,                // solver module
                        OpenWQ_output,
                        openwq_source_EWF_name,
                        irow, icol,
                        inflowi);

                }
            }     
        }
    }
    else{
        std::cout << "inflow ix,iy in NODATA_VALUE cell" << std::endl; // err message
        errflag = true;
    }
    
}catch (int e){

    std::cout << "problem in 'add_inflow' module" << std::endl; // err message
    errflag = true;
}

return errflag;

}

// Check if directory exists and create it
void check_mkdir(
    std::string &dirname){
    
    struct stat st = {0};
    
    // convert to *char
    const char *cstr = dirname.c_str();

    // mkdir
    if (stat(cstr, &st) == -1) {
        mkdir(cstr, 0700);
    }

}