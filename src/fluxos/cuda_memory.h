// Copyright 2019, Diogo Costa
// CUDA GPU acceleration for FLUXOS

#ifndef CUDA_MEMORY_H_INCLUDED
#define CUDA_MEMORY_H_INCLUDED

#ifdef USE_CUDA

#include <cstddef>

class GlobVar;

struct CudaGridData {
    size_t MROWS, MCOLS;
    size_t total_size; // MROWS * MCOLS

    // Scalar constants (updated each timestep)
    double hdry, gacc, cfl, cvdef, nuem, dtfl;
    size_t dxy, arbase;
    unsigned int NROWS, NCOLS;
    double NODATA_VALUE;

    // Device pointers for double matrices
    double *d_z, *d_zb, *d_h, *d_ux, *d_uy;
    double *d_qx, *d_qy, *d_qxf, *d_qyf;
    double *d_us, *d_dh, *d_dqx, *d_dqy;
    double *d_ks;
    double *d_fe_1, *d_fe_2, *d_fe_3;
    double *d_fn_1, *d_fn_2, *d_fn_3;
    double *d_twetimetracer;
    double *d_h0;

    // Device pointers for float matrices
    float *d_ldry, *d_innerNeumannBCWeir, *d_ldry_prev;

    // Device pointer for concentration (one species at a time)
    double *d_conc_SW;

    // Soil infiltration device arrays (Horton decay model)
    double *d_soil_infil_rate;
    double *d_soil_Ks;
    double *d_soil_f0;
    double *d_soil_k;
    double *d_soil_wetting_time;
    bool soil_allocated;

    // Reduction buffers
    double *d_block_reduce;
    double *d_reduce_result;   // preallocated 2-element buffer for final reduction
    size_t num_blocks;
};

class CudaMemoryManager {
public:
    CudaGridData grid;

    CudaMemoryManager();
    ~CudaMemoryManager();

    void allocate(size_t MROWS, size_t MCOLS);
    void deallocate();

    // Bulk transfer: all fields from Armadillo to GPU (call once at init)
    void copy_all_to_device(GlobVar& ds);

    // Bulk transfer: all fields from GPU back to Armadillo
    void copy_all_to_host(GlobVar& ds);

    // Transfer only output-relevant fields (for write_results)
    void copy_output_fields_to_host(GlobVar& ds);

    // Update scalar constants on device (call each timestep)
    void update_scalars(GlobVar& ds);

    // Concentration field for a specific chemical species
    void copy_conc_to_device(GlobVar& ds, int ichem);
    void copy_conc_to_host(GlobVar& ds, int ichem);

    // Soil infiltration GPU memory
    void allocate_soil();
    void copy_soil_to_device(GlobVar& ds);
    void copy_soil_to_host(GlobVar& ds);

    void sync();

private:
    bool allocated;
};

#endif // USE_CUDA
#endif // CUDA_MEMORY_H_INCLUDED
