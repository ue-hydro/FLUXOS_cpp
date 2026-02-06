cmake_minimum_required(VERSION 3.18)

# ============================================================================
# FLUXOS-OVERLAND CMake Configuration
# Hybrid MPI+OpenMP parallelization for HPC clusters
# Optional CUDA GPU acceleration
# ============================================================================
#
# Build Options:
#   -DMODE_release=ON     Enable release mode with optimizations (default: OFF)
#   -DUSE_MPI=ON          Enable MPI for distributed computing (default: OFF)
#   -DUSE_CUDA=ON         Enable CUDA GPU acceleration (default: OFF)
#   -DUSE_TRIMESH=ON      Enable triangular mesh support (default: OFF)
#   -DUSE_OPENWQ=ON       Enable OpenWQ water quality coupling (default: OFF)
#
# Build Examples:
#   cmake -DMODE_release=ON ..                              # OpenMP only (single node)
#   cmake -DMODE_release=ON -DUSE_MPI=ON ..                 # MPI+OpenMP (multi-node HPC)
#   cmake -DMODE_release=ON -DUSE_CUDA=ON ..                # OpenMP+CUDA (single GPU)
#   cmake -DMODE_release=ON -DUSE_TRIMESH=ON ..             # OpenMP + triangular mesh
#   cmake -DMODE_release=ON -DUSE_TRIMESH=ON -DUSE_CUDA=ON ..  # Triangular + CUDA
#
# ============================================================================
# Scalability Guidelines
# ============================================================================
#
#   Domain Size         Recommended Configuration
#   ---------------     --------------------------------------------------
#   < 1000 x 1000       OpenMP only (single node)
#   1000 - 5000²        4-16 MPI processes
#   5000 - 10000²       16-64 MPI processes
#   > 10000²            64+ MPI processes
#
# For hybrid MPI+OpenMP, use 2-4 OpenMP threads per MPI process.
# Example for 64 cores: 16 MPI processes x 4 OpenMP threads
#
# ============================================================================

# Define compilation mode
option(MODE_release "ON: release mode; OFF: debug mode" OFF)

# Define MPI option for HPC clusters
option(USE_MPI "ON: Enable MPI for distributed computing; OFF: Single node only" OFF)

# Define CUDA option for GPU acceleration
option(USE_CUDA "ON: Enable CUDA GPU acceleration; OFF: CPU only" OFF)

# Define triangular mesh option
option(USE_TRIMESH "ON: Enable unstructured triangular mesh support; OFF: Regular mesh only" OFF)

# Define OpenWQ water quality coupling option
option(USE_OPENWQ "ON: Enable OpenWQ water quality coupling; OFF: Standalone hydro" OFF)

# Set compilation mode flags
IF(MODE_release MATCHES ON)
        MESSAGE(STATUS "compilation model : release")
        # Aggressive optimization flags for maximum performance
        # Note: Using -fno-finite-math-only to maintain IEEE compliance for JSON library
        set(CMAKE_CXX_FLAGS "-O3 -march=native -mtune=native -funroll-loops -ftree-vectorize -fno-math-errno")
        # Enable link-time optimization for cross-file inlining
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
        set(FLUXOS_EXEC_name "fluxos")
ELSE(MODE_release MATCHES OFF)
        MESSAGE(STATUS "compilation model : debug")
        set(CMAKE_CXX_FLAGS "-g3 -Wall -pedantic")
        set(FLUXOS_EXEC_name "fluxos_debug")
ENDIF()

# Settings depending on OS
IF (UNIX)
	IF (UNIX AND NOT APPLE)
		MESSAGE(STATUS "current platform : Linux")
	ELSEIF(APPLE)
		MESSAGE(STATUS "current platform : MacOS")
	ENDIF()

	project(fluxos)

	file(GLOB SOURCES "src/fluxos/*.cpp")

	# OpenWQ water quality coupling
	# Always compile the OpenWQ hydrolink bridge file.
	# When USE_OPENWQ=OFF, stub implementations are used.
	# When USE_OPENWQ=ON, the full OpenWQ library must be available.
	file(GLOB OPENWQ_SOURCES "src/openwq/*.cpp")
	IF(USE_OPENWQ MATCHES ON)
		MESSAGE(STATUS "OpenWQ support: ENABLED")
		add_definitions(-DUSE_OPENWQ)
	ELSE()
		MESSAGE(STATUS "OpenWQ support: DISABLED (use -DUSE_OPENWQ=ON to enable)")
	ENDIF()

	set(CMAKE_CXX_STANDARD 17)

	# Packages
	find_package(Armadillo REQUIRED)
	FIND_PACKAGE(HDF5 REQUIRED)
	find_package(OpenMP)

	# MPI support for HPC clusters
	IF(USE_MPI MATCHES ON)
		find_package(MPI REQUIRED)
		MESSAGE(STATUS "MPI support : ENABLED")
		MESSAGE(STATUS "MPI C compiler: ${MPI_C_COMPILER}")
		MESSAGE(STATUS "MPI CXX compiler: ${MPI_CXX_COMPILER}")
		add_definitions(-DUSE_MPI)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
		# Update executable name to indicate MPI version
		IF(MODE_release MATCHES ON)
			set(FLUXOS_EXEC_name "fluxos_mpi")
		ELSE()
			set(FLUXOS_EXEC_name "fluxos_mpi_debug")
		ENDIF()
	ELSE()
		MESSAGE(STATUS "MPI support : DISABLED (use -DUSE_MPI=ON to enable)")
	ENDIF()

	# CUDA GPU acceleration support
	IF(USE_CUDA MATCHES ON)
		enable_language(CUDA)
		find_package(CUDAToolkit REQUIRED)
		MESSAGE(STATUS "CUDA support : ENABLED")
		MESSAGE(STATUS "CUDA compiler: ${CMAKE_CUDA_COMPILER}")
		MESSAGE(STATUS "CUDA toolkit : ${CUDAToolkit_VERSION}")
		add_definitions(-DUSE_CUDA)

		# Set CUDA architectures (user can override with -DCMAKE_CUDA_ARCHITECTURES=xx)
		if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
			set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86 89 90)
		endif()

		set(CMAKE_CUDA_STANDARD 17)

		# CUDA compile flags
		IF(MODE_release MATCHES ON)
			set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --use_fast_math --expt-relaxed-constexpr")
		ELSE()
			set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -G --expt-relaxed-constexpr")
		ENDIF()

		# Collect CUDA source files
		file(GLOB CUDA_SOURCES "src/fluxos/*.cu")

		# Update executable name
		IF(MODE_release MATCHES ON)
			set(FLUXOS_EXEC_name "fluxos_cuda")
		ELSE()
			set(FLUXOS_EXEC_name "fluxos_cuda_debug")
		ENDIF()
	ELSE()
		MESSAGE(STATUS "CUDA support : DISABLED (use -DUSE_CUDA=ON to enable)")
		set(CUDA_SOURCES "")
	ENDIF()

	# Triangular mesh support
	IF(USE_TRIMESH MATCHES ON)
		MESSAGE(STATUS "TriMesh support: ENABLED")
		add_definitions(-DUSE_TRIMESH)

		# Collect triangular mesh C++ sources
		file(GLOB TRIMESH_SOURCES "src/fluxos/trimesh/*.cpp")

		# Collect triangular mesh CUDA sources if CUDA is enabled
		IF(USE_CUDA MATCHES ON)
			file(GLOB TRIMESH_CUDA_SOURCES "src/fluxos/trimesh/*.cu")
		ELSE()
			set(TRIMESH_CUDA_SOURCES "")
		ENDIF()
	ELSE()
		MESSAGE(STATUS "TriMesh support: DISABLED (use -DUSE_TRIMESH=ON to enable)")
		set(TRIMESH_SOURCES "")
		set(TRIMESH_CUDA_SOURCES "")
	ENDIF()

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

	include_directories(${AMRADILLO_INCLUDE_DIR} ${PROJECT_SOURCE_DIR})

	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

	# Create Results folder if non-existant
	file(MAKE_DIRECTORY Results)

	add_executable(${FLUXOS_EXEC_name}
			${SOURCES}
			${OPENWQ_SOURCES}
			${CUDA_SOURCES}
			${TRIMESH_SOURCES}
			${TRIMESH_CUDA_SOURCES})

	target_include_directories(
		${FLUXOS_EXEC_name} PUBLIC
		${HDF5_INCLUDE_DIRS}
		${ARMADILLO_INCLUDE_DIR})

	# Add MPI include directories if enabled
	IF(USE_MPI MATCHES ON)
		target_include_directories(
			${FLUXOS_EXEC_name} PUBLIC
			${MPI_CXX_INCLUDE_DIRS})
	ENDIF()

	target_link_libraries(
		${FLUXOS_EXEC_name}
		${HDF5_C_LIBRARY_hdf5}
		${ARMADILLO_LIBRARIES}
		OpenMP::OpenMP_CXX)

	# Add MPI libraries if enabled
	IF(USE_MPI MATCHES ON)
		target_link_libraries(
			${FLUXOS_EXEC_name}
			${MPI_CXX_LIBRARIES})
	ENDIF()

	# Add CUDA libraries if enabled
	IF(USE_CUDA MATCHES ON)
		target_link_libraries(
			${FLUXOS_EXEC_name}
			CUDA::cudart)
		set_target_properties(${FLUXOS_EXEC_name} PROPERTIES
			CUDA_SEPARABLE_COMPILATION ON)
	ENDIF()


ELSE()
	MESSAGE(STATUS "curent platform : windows")
	project(fluxos)
	file(GLOB SOURCES "src/*.cpp")
	# Create Results folder if non-existant
	file(MAKE_DIRECTORY Results)

	set(CMAKE_CXX_STANDARD 17)

	# Packages
	find_package(OpenMP)

	# MPI support for Windows HPC
	IF(USE_MPI MATCHES ON)
		find_package(MPI REQUIRED)
		MESSAGE(STATUS "MPI support : ENABLED")
		add_definitions(-DUSE_MPI)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
	ELSE()
		MESSAGE(STATUS "MPI support : DISABLED")
	ENDIF()

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Results)
	add_executable(${FLUXOS_EXEC_name}, ${SOURCES})
	target_include_directories(
		${FLUXOS_EXEC_name} PUBLIC
		${HDF5_INCLUDE_DIRS}
		"/usr/include/")

	IF(USE_MPI MATCHES ON)
		target_include_directories(
			${FLUXOS_EXEC_name} PUBLIC
			${MPI_CXX_INCLUDE_DIRS})
		target_link_libraries(
			${FLUXOS_EXEC_name}
			${MPI_CXX_LIBRARIES})
	ENDIF()

ENDIF()
