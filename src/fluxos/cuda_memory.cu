// Copyright 2019, Diogo Costa
// CUDA GPU acceleration for FLUXOS - Memory management

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <armadillo>
#include "GlobVar.h"
#include "cuda_memory.h"

#define CUDA_CHECK(call)                                                       \
    do {                                                                        \
        cudaError_t err = call;                                                 \
        if (err != cudaSuccess) {                                               \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,   \
                    cudaGetErrorString(err));                                    \
            exit(EXIT_FAILURE);                                                 \
        }                                                                       \
    } while (0)

// Helper: allocate a device double array
static void alloc_d(double** ptr, size_t n) {
    CUDA_CHECK(cudaMalloc(ptr, n * sizeof(double)));
    CUDA_CHECK(cudaMemset(*ptr, 0, n * sizeof(double)));
}

// Helper: allocate a device float array
static void alloc_f(float** ptr, size_t n) {
    CUDA_CHECK(cudaMalloc(ptr, n * sizeof(float)));
    CUDA_CHECK(cudaMemset(*ptr, 0, n * sizeof(float)));
}

// Helper: copy Armadillo double mat to device (column-major, direct memcpy)
static void arma_to_device(const arma::Mat<double>& mat, double* d_ptr, size_t n) {
    CUDA_CHECK(cudaMemcpy(d_ptr, mat.memptr(), n * sizeof(double), cudaMemcpyHostToDevice));
}

// Helper: copy Armadillo float mat to device
static void arma_to_device(const arma::Mat<float>& mat, float* d_ptr, size_t n) {
    CUDA_CHECK(cudaMemcpy(d_ptr, mat.memptr(), n * sizeof(float), cudaMemcpyHostToDevice));
}

// Helper: copy device double array back to Armadillo mat
static void device_to_arma(const double* d_ptr, arma::Mat<double>& mat, size_t n) {
    CUDA_CHECK(cudaMemcpy(mat.memptr(), d_ptr, n * sizeof(double), cudaMemcpyDeviceToHost));
}

// Helper: copy device float array back to Armadillo mat
static void device_to_arma(const float* d_ptr, arma::Mat<float>& mat, size_t n) {
    CUDA_CHECK(cudaMemcpy(mat.memptr(), d_ptr, n * sizeof(float), cudaMemcpyDeviceToHost));
}

CudaMemoryManager::CudaMemoryManager() : allocated(false) {
    memset(&grid, 0, sizeof(CudaGridData));
}

CudaMemoryManager::~CudaMemoryManager() {
    if (allocated) deallocate();
}

void CudaMemoryManager::allocate(size_t MROWS, size_t MCOLS) {
    grid.MROWS = MROWS;
    grid.MCOLS = MCOLS;
    grid.total_size = MROWS * MCOLS;
    size_t n = grid.total_size;

    // Print GPU info
    int device;
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDevice(&device));
    CUDA_CHECK(cudaGetDeviceProperties(&prop, device));
    printf("CUDA device: %s (%.1f GB, compute %d.%d)\n",
           prop.name, prop.totalGlobalMem / 1e9, prop.major, prop.minor);

    size_t total_bytes = 0;

    // Double matrices (23)
    alloc_d(&grid.d_z, n);   total_bytes += n * sizeof(double);
    alloc_d(&grid.d_zb, n);  total_bytes += n * sizeof(double);
    alloc_d(&grid.d_h, n);   total_bytes += n * sizeof(double);
    alloc_d(&grid.d_ux, n);  total_bytes += n * sizeof(double);
    alloc_d(&grid.d_uy, n);  total_bytes += n * sizeof(double);
    alloc_d(&grid.d_qx, n);  total_bytes += n * sizeof(double);
    alloc_d(&grid.d_qy, n);  total_bytes += n * sizeof(double);
    alloc_d(&grid.d_qxf, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_qyf, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_us, n);  total_bytes += n * sizeof(double);
    alloc_d(&grid.d_dh, n);  total_bytes += n * sizeof(double);
    alloc_d(&grid.d_dqx, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_dqy, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_ks, n);  total_bytes += n * sizeof(double);
    alloc_d(&grid.d_fe_1, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_fe_2, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_fe_3, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_fn_1, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_fn_2, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_fn_3, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_twetimetracer, n); total_bytes += n * sizeof(double);
    alloc_d(&grid.d_h0, n);         total_bytes += n * sizeof(double);

    // Float matrices (3)
    alloc_f(&grid.d_ldry, n);              total_bytes += n * sizeof(float);
    alloc_f(&grid.d_innerNeumannBCWeir, n); total_bytes += n * sizeof(float);
    alloc_f(&grid.d_ldry_prev, n);          total_bytes += n * sizeof(float);

    // Concentration (one species buffer)
    alloc_d(&grid.d_conc_SW, n); total_bytes += n * sizeof(double);

    // Reduction buffer: max blocks for a 16x16 launch
    unsigned int NROWS_est = MROWS > 2 ? MROWS - 2 : 1;
    unsigned int NCOLS_est = MCOLS > 2 ? MCOLS - 2 : 1;
    grid.num_blocks = ((NROWS_est + 15) / 16) * ((NCOLS_est + 15) / 16);
    alloc_d(&grid.d_block_reduce, grid.num_blocks * 2); // *2 for min+max
    total_bytes += grid.num_blocks * 2 * sizeof(double);

    printf("CUDA memory allocated: %.1f MB for %zux%zu grid\n",
           total_bytes / 1e6, MROWS, MCOLS);

    allocated = true;
}

void CudaMemoryManager::deallocate() {
    if (!allocated) return;

    CUDA_CHECK(cudaFree(grid.d_z));
    CUDA_CHECK(cudaFree(grid.d_zb));
    CUDA_CHECK(cudaFree(grid.d_h));
    CUDA_CHECK(cudaFree(grid.d_ux));
    CUDA_CHECK(cudaFree(grid.d_uy));
    CUDA_CHECK(cudaFree(grid.d_qx));
    CUDA_CHECK(cudaFree(grid.d_qy));
    CUDA_CHECK(cudaFree(grid.d_qxf));
    CUDA_CHECK(cudaFree(grid.d_qyf));
    CUDA_CHECK(cudaFree(grid.d_us));
    CUDA_CHECK(cudaFree(grid.d_dh));
    CUDA_CHECK(cudaFree(grid.d_dqx));
    CUDA_CHECK(cudaFree(grid.d_dqy));
    CUDA_CHECK(cudaFree(grid.d_ks));
    CUDA_CHECK(cudaFree(grid.d_fe_1));
    CUDA_CHECK(cudaFree(grid.d_fe_2));
    CUDA_CHECK(cudaFree(grid.d_fe_3));
    CUDA_CHECK(cudaFree(grid.d_fn_1));
    CUDA_CHECK(cudaFree(grid.d_fn_2));
    CUDA_CHECK(cudaFree(grid.d_fn_3));
    CUDA_CHECK(cudaFree(grid.d_twetimetracer));
    CUDA_CHECK(cudaFree(grid.d_h0));
    CUDA_CHECK(cudaFree(grid.d_ldry));
    CUDA_CHECK(cudaFree(grid.d_innerNeumannBCWeir));
    CUDA_CHECK(cudaFree(grid.d_ldry_prev));
    CUDA_CHECK(cudaFree(grid.d_conc_SW));
    CUDA_CHECK(cudaFree(grid.d_block_reduce));

    allocated = false;
    printf("CUDA memory freed\n");
}

void CudaMemoryManager::copy_all_to_device(GlobVar& ds) {
    size_t n = grid.total_size;

    arma_to_device(*ds.z, grid.d_z, n);
    arma_to_device(*ds.zb, grid.d_zb, n);
    arma_to_device(*ds.h, grid.d_h, n);
    arma_to_device(*ds.ux, grid.d_ux, n);
    arma_to_device(*ds.uy, grid.d_uy, n);
    arma_to_device(*ds.qx, grid.d_qx, n);
    arma_to_device(*ds.qy, grid.d_qy, n);
    arma_to_device(*ds.qxf, grid.d_qxf, n);
    arma_to_device(*ds.qyf, grid.d_qyf, n);
    arma_to_device(*ds.us, grid.d_us, n);
    arma_to_device(*ds.dh, grid.d_dh, n);
    arma_to_device(*ds.dqx, grid.d_dqx, n);
    arma_to_device(*ds.dqy, grid.d_dqy, n);
    arma_to_device(*ds.ks, grid.d_ks, n);
    arma_to_device(*ds.fe_1, grid.d_fe_1, n);
    arma_to_device(*ds.fe_2, grid.d_fe_2, n);
    arma_to_device(*ds.fe_3, grid.d_fe_3, n);
    arma_to_device(*ds.fn_1, grid.d_fn_1, n);
    arma_to_device(*ds.fn_2, grid.d_fn_2, n);
    arma_to_device(*ds.fn_3, grid.d_fn_3, n);
    arma_to_device(*ds.twetimetracer, grid.d_twetimetracer, n);
    arma_to_device(*ds.h0, grid.d_h0, n);
    arma_to_device(*ds.ldry, grid.d_ldry, n);
    arma_to_device(*ds.innerNeumannBCWeir, grid.d_innerNeumannBCWeir, n);
    arma_to_device(*ds.ldry_prev, grid.d_ldry_prev, n);

    // Copy scalars
    update_scalars(ds);

    printf("CUDA: all fields copied to device\n");
}

void CudaMemoryManager::copy_all_to_host(GlobVar& ds) {
    size_t n = grid.total_size;

    device_to_arma(grid.d_z, *ds.z, n);
    device_to_arma(grid.d_zb, *ds.zb, n);
    device_to_arma(grid.d_h, *ds.h, n);
    device_to_arma(grid.d_ux, *ds.ux, n);
    device_to_arma(grid.d_uy, *ds.uy, n);
    device_to_arma(grid.d_qx, *ds.qx, n);
    device_to_arma(grid.d_qy, *ds.qy, n);
    device_to_arma(grid.d_qxf, *ds.qxf, n);
    device_to_arma(grid.d_qyf, *ds.qyf, n);
    device_to_arma(grid.d_us, *ds.us, n);
    device_to_arma(grid.d_dh, *ds.dh, n);
    device_to_arma(grid.d_dqx, *ds.dqx, n);
    device_to_arma(grid.d_dqy, *ds.dqy, n);
    device_to_arma(grid.d_ks, *ds.ks, n);
    device_to_arma(grid.d_fe_1, *ds.fe_1, n);
    device_to_arma(grid.d_fe_2, *ds.fe_2, n);
    device_to_arma(grid.d_fe_3, *ds.fe_3, n);
    device_to_arma(grid.d_fn_1, *ds.fn_1, n);
    device_to_arma(grid.d_fn_2, *ds.fn_2, n);
    device_to_arma(grid.d_fn_3, *ds.fn_3, n);
    device_to_arma(grid.d_twetimetracer, *ds.twetimetracer, n);
    device_to_arma(grid.d_h0, *ds.h0, n);
    device_to_arma(grid.d_ldry, *ds.ldry, n);
    device_to_arma(grid.d_innerNeumannBCWeir, *ds.innerNeumannBCWeir, n);
    device_to_arma(grid.d_ldry_prev, *ds.ldry_prev, n);
}

void CudaMemoryManager::copy_output_fields_to_host(GlobVar& ds) {
    size_t n = grid.total_size;

    device_to_arma(grid.d_z, *ds.z, n);
    device_to_arma(grid.d_h, *ds.h, n);
    device_to_arma(grid.d_qx, *ds.qx, n);
    device_to_arma(grid.d_qy, *ds.qy, n);
    device_to_arma(grid.d_us, *ds.us, n);
    device_to_arma(grid.d_ldry, *ds.ldry, n);
    device_to_arma(grid.d_twetimetracer, *ds.twetimetracer, n);
    device_to_arma(grid.d_h0, *ds.h0, n);
    device_to_arma(grid.d_ldry_prev, *ds.ldry_prev, n);
    device_to_arma(grid.d_fe_1, *ds.fe_1, n);
    device_to_arma(grid.d_fn_1, *ds.fn_1, n);
}

void CudaMemoryManager::update_scalars(GlobVar& ds) {
    grid.NROWS = ds.NROWS;
    grid.NCOLS = ds.NCOLS;
    grid.hdry = ds.hdry;
    grid.gacc = ds.gacc;
    grid.cfl = ds.cfl;
    grid.cvdef = ds.cvdef;
    grid.nuem = ds.nuem;
    grid.dtfl = ds.dtfl;
    grid.dxy = ds.dxy;
    grid.arbase = ds.arbase;
    grid.NODATA_VALUE = ds.NODATA_VALUE;
}

void CudaMemoryManager::copy_conc_to_device(GlobVar& ds, int ichem) {
    arma::Mat<double>& conc = (*ds.conc_SW)[ichem];
    CUDA_CHECK(cudaMemcpy(grid.d_conc_SW, conc.memptr(),
                          grid.total_size * sizeof(double), cudaMemcpyHostToDevice));
}

void CudaMemoryManager::copy_conc_to_host(GlobVar& ds, int ichem) {
    arma::Mat<double>& conc = (*ds.conc_SW)[ichem];
    CUDA_CHECK(cudaMemcpy(conc.memptr(), grid.d_conc_SW,
                          grid.total_size * sizeof(double), cudaMemcpyDeviceToHost));
}

void CudaMemoryManager::sync() {
    CUDA_CHECK(cudaDeviceSynchronize());
}

#endif // USE_CUDA
