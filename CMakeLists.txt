
#  CMakeLists.txt
#  
#  NiftyRec
#  Stefano Pedemonte, May 2012.
#  CMIC - Centre for Medical Image Computing 
#  UCL - University College London. 
#  Released under BSD licence, see LICENSE.txt 

#-----------------------------------------------------------------------------
PROJECT(nifty_rec)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
cmake_policy(VERSION 2.8)
INCLUDE(ExternalProject)
ENABLE_TESTING()

#-----------------------------------------------------------------------------
SET(Nifty_Rec_VERSION_MAJOR 1)
SET(Nifty_Rec_VERSION_MINOR 6)

#-----------------------------------------------------------------------------
OPTION(NiftyRec_BUILD_APPS          "Build command line applications (including 2D Slicer plugins)" ON)
OPTION(NiftyRec_BUILD_MEX           "Build Matlab extension modules" ON)
OPTION(NiftyRec_BUILD_PY            "Build Python extension modules" ON)
OPTION(NiftyRec_USE_CUDA            "Compile GPU implementation of the code" OFF)
IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    IF(NOT(${CMAKE_SYSTEM_NAME} MATCHES "Windows"))
        IF(NOT CMAKE_SYSTEM_VERSION VERSION_GREATER 11.0 OR CMAKE_SYSTEM_VERSION VERSION_EQUAL 11.0.0)
            OPTION(NiftyRec_USE_OPENMP "To use openMP for multi-CPU processing" ON)
        ENDIF()
    ENDIF()
ENDIF()
OPTION(NiftyRec_SUPPORT_NRRD        "Support NRRD image file format. Enable to generate 3D Slicer plugins." OFF)
OPTION(NiftyRec_COMPILE_DOC         "Compile pdf documentation from the Latex source" OFF)
OPTION(NiftyRec_INCLUDE_DATA        "Include data sets for examples and demos" OFF)
OPTION(NiftyRec_BUILD_RENDERER      "Build multi-modal volume renderer" OFF)
OPTION(NiftyRec_VERBOSE             "Print out extra information when compiling" OFF)
#OPTION(NiftyRec_INSTALL_NIFTYSEG    "Download and install NiftySeg" OFF)
#OPTION(NiftyRec_BUILD_SHARED        "Build as shared libraries" ON)

IF(NOT(${CMAKE_SYSTEM_NAME} MATCHES "Windows"))
    SET(NiftyRec_BUILD_SHARED ON)
ELSE()
    SET(NiftyRec_BUILD_SHARED OFF)
ENDIF()
SET(NiftyRec_DO_PACKAGE ON)

IF(CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs" FORCE)
ENDIF(CMAKE_CONFIGURATION_TYPES)
SET(CMAKE_BUILD_TYPE "Release")

#-----------------------------------------------------------------------------

IF(NiftyRec_BUILD_PY)
    FIND_PACKAGE(PythonInterp)
ENDIF(NiftyRec_BUILD_PY)
#Find NiftyReg (commented out as NiftyReg libraries have been temporarily included in the source tree)
#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/nifty_reg")
#FIND_PACKAGE(NIFTY_REG REQUIRED)
#INCLUDE_DIRECTORIES(${NIFTY_REG_INCLUDE_DIRS}) #Nifty-Reg libraries are in NIFTY_REG_LIBRARY

#-----------------------------------------------------------------------------

SET(LIB_INSTALL_DIR "lib")
SET(BIN_INSTALL_DIR "niftyrec/bin")
SET(INC_INSTALL_DIR "include")
SET(DOC_INSTALL_DIR "niftyrec/documentation")
SET(MEX_INSTALL_DIR_NIFTYREC "niftyrec/matlab")
SET(MEX_INSTALL_DIR_NIFTYREG "niftyreg/matlab")
SET(MEX_INSTALL_DIR_NIFTYSEG "niftyseg/matlab")
SET(PYTHON_INSTALL_DIR "niftyrec/python")

#-----------------------------------------------------------------------------

#Point rpath of all libraries to lib install dir 
#This avoids setting LD_LIBRARY_PATH 
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR})

#-----------------------------------------------------------------------------
FIND_PACKAGE(Threads)

IF(NiftyRec_USE_OPENMP)
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        IF(CMAKE_SYSTEM_VERSION VERSION_GREATER 11.0 OR CMAKE_SYSTEM_VERSION VERSION_EQUAL 11.0.0)
            MESSAGE("OpenMP appears to be broken using Mac OS Lion")
            MESSAGE(FATAL_ERROR "Please turn NiftyRec_USE_OPENMP OFF")
        ENDIF()
    ENDIF()
	
    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)

    IF(GCC_VERSION VERSION_GREATER 4.2 OR GCC_VERSION VERSION_EQUAL 4.2)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    ENDIF()
	
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
    ENDIF()
ENDIF(NiftyRec_USE_OPENMP)

#-----------------------------------------------------------------------------

IF(UNIX AND NiftyRec_BUILD_MEX AND NOT NiftyRec_BUILD_SHARED)
        MESSAGE("Unix platform: desabling Matlab extension because libraries are not built as shared.")
        SET(NiftyRec_BUILD_MEX OFF)
ENDIF()

IF(UNIX AND NiftyRec_BUILD_PY AND NOT NiftyRec_BUILD_SHARED)
        MESSAGE("Unix platform: desabling Python extension because libraries are not built as shared.")
        SET(NiftyRec_BUILD_PY OFF)
ENDIF()
IF(NiftyRec_BUILD_PY)
        IF(NOT PYTHONINTERP_FOUND)
                MESSAGE("Python interpreter not found")
                SET(NiftyRec_BUILD_PY OFF)        
        ENDIF()
ENDIF(NiftyRec_BUILD_PY)

#-----------------------------------------------------------------------------

IF(NiftyRec_BUILD_SHARED)
        SET(LIB_TYPE SHARED)
        IF(UNIX)
                IF("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
                        ADD_DEFINITIONS(-fPIC)
                        SET(CMAKE_CC_FLAGS_RELEASE "${CMAKE_CC_FLAGS_RELEASE} -fPIC -ldl"  )
                        SET(CMAKE_CC_FLAGS_DEBUG "${CMAKE_CC_FLAGS_DEBUG} -fPIC -ldl"  )
                        SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC -ldl"  )
                        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fPIC -ldl"  )
                ENDIF("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
        ENDIF(UNIX)
ELSE(NiftyRec_BUILD_SHARED)
        SET(LIB_TYPE STATIC)
ENDIF(NiftyRec_BUILD_SHARED)

#-----------------------------------------------------------------------------

IF(NiftyRec_VERBOSE)
	ADD_DEFINITIONS(-D_VERBOSE)
ENDIF(NiftyRec_VERBOSE)

SET(USE_CUDA 0)
IF(NiftyRec_USE_CUDA)
        SET(USE_CUDA 1)
        FIND_PACKAGE(CUDA REQUIRED)
        ADD_DEFINITIONS(-D_USE_CUDA)
        SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-D_USE_CUDA")
        SET(CUDA_CUT_INCLUDE "${CUDA_SDK_ROOT_DIR}/common/inc/")
        INCLUDE_DIRECTORIES(${CUDA_CUT_INCLUDE})
        INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})
        MESSAGE("Cuda toolkit root:    " ${CUDA_TOOLKIT_ROOT_DIR})
        MESSAGE("Cuda sdk root:        " ${CUDA_SDK_ROOT_DIR})
        MESSAGE("Cuda toolkit include: " ${CUDA_TOOLKIT_INCLUDE})
        MESSAGE("Cuda cut include:     " ${CUDA_CUT_INCLUDE})
        CUDA_INCLUDE_DIRECTORIES(${CUDA_CUT_INCLUDE_DIR})
ENDIF(NiftyRec_USE_CUDA)

IF(NiftyRec_SUPPORT_NRRD)
        ADD_DEFINITIONS(-D_SUPPORT_NRRD)
ENDIF(NiftyRec_SUPPORT_NRRD)

#-----------------------------------------------------------------------------

#Build the following libraries. Names are set here (not in the subfolders) so targets 
#in the subfolders know what to link against. In order to add a new library, list it here. 
SET(NIFTY_REC_ET_INTERF_LIBRARY _et_array_interface)
SET(NIFTY_REC_TT_INTERF_LIBRARY _tt_array_interface)
SET(NIFTY_REC_REG_INTERF_LIBRARY _reg_array_interface)
SET(NIFTY_REC_ET_LIBRARIES _et_common _niftyrec_memory _et_line_integral _et_line_integral_attenuated _et_line_backproject_attenuated _et_line_backproject _et_accumulate _et_clear_accumulator _et_convolve2D _et_convolveSeparable2D _et_convolveFFT2D)
SET(NIFTY_REC_ET_LIBRARIES_GPU _et_line_integral_gpu _et_line_integral_attenuated_gpu _et_line_backproject_gpu _et_line_backproject_attenuated_gpu _et_clear_accumulator_gpu _et_accumulate_gpu _et_convolveSeparable2D_gpu _et_convolveFFT2D_gpu _et_attenuation_gradient_gpu)
SET(NIFTY_REC_ET_MEX et_project_mex et_project_partial_mex et_backproject_mex et_project_backproject_mex et_convolve_mex et_rotate_mex et_list_gpus_mex et_set_gpu_mex et_affine_mex et_reorder_activity_in et_reorder_activity_out et_variance_centre_approx_mex et_fisher_grid_mex et_fisher_grid_invprojection_mex et_gradient_attenuation_mex et_reset_gpu_mex)
SET(NIFTY_REC_TT_LIBRARIES _tt_common _tt_line_project_ray_cpu _tt_line_backproject_ray_cpu)
SET(NIFTY_REC_TT_LIBRARIES_GPU _tt_line_project_ray_gpu _tt_line_backproject_ray_gpu) 
SET(NIFTY_REC_TT_MEX tt_project_ray_mex tt_backproject_ray_mex) 

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/emission-lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/transmission-lib)



IF(NiftyRec_USE_CUDA)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/emission-lib_gpu)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/transmission-lib_gpu)
ENDIF(NiftyRec_USE_CUDA)

#Build NiftyReg Libraries (they have been added to the source tree to remove 
#dependency - there is currently no binary installer for NiftyReg)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/nifti)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/reg-lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/reg-lib_gpu)
SUBDIRS(nifti)
SUBDIRS(tclap)
IF(NiftyRec_SUPPORT_NRRD)
        SUBDIRS(teem)
ENDIF(NiftyRec_SUPPORT_NRRD)
SUBDIRS(reg-lib)
SUBDIRS(reg-lib_gpu)
SET(NIFTY_REG_LIBRARIES reg_nifti _reg_tools _reg_resampling _reg_affineTransformation _reg_bspline _reg_ssd _reg_mutualinformation _reg_blockMatching)
SET(NIFTY_REG_LIBRARIES_GPU _reg_cudaCommon _reg_resampling_gpu _reg_affineTransformation_gpu _reg_mutualinformation_gpu _reg_bspline_gpu _reg_tools_gpu _reg_blockMatching_gpu)


SUBDIRS(emission-lib)
SUBDIRS(transmission-lib)
IF(NiftyRec_BUILD_MEX)
        SUBDIRS(Matlab)
ENDIF(NiftyRec_BUILD_MEX)
IF(NiftyRec_BUILD_PY)
        SUBDIRS(Python)
ENDIF(NiftyRec_BUILD_PY)
IF(NiftyRec_BUILD_RENDERER)
        SUBDIRS(volume_render)
ENDIF(NiftyRec_BUILD_RENDERER)
IF(NiftyRec_USE_CUDA)
	SUBDIRS(emission-lib_gpu)
	SUBDIRS(transmission-lib_gpu)
ENDIF(NiftyRec_USE_CUDA)
IF(NiftyRec_COMPILE_DOC)
        SUBDIRS(documentation)
ENDIF(NiftyRec_COMPILE_DOC)
IF(NiftyRec_BUILD_APPS)
        SUBDIRS(apps)
ENDIF(NiftyRec_BUILD_APPS)
SUBDIRS(misc)
IF(NiftyRec_INCLUDE_DATA)
    SUBDIRS(data)
ENDIF(NiftyRec_INCLUDE_DATA)

#-----------------------------------------------------------------------------
#IF(NiftyRec_INSTALL_NIFTYSEG)
#ExternalProject_Add(niftyseg
#            SOURCE_DIR ${PROJECT_BINARY_DIR}/niftyseg
#            PREFIX ${PROJECT_BINARY_DIR}/niftyseg
#            STAMP_DIR ${PROJECT_BINARY_DIR}/niftyseg/stamps
#            TMP_DIR ${PROJECT_BINARY_DIR}/niftyseg/stamps
#            DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/niftyseg
#            BINARY_DIR ${PROJECT_BINARY_DIR}/niftyseg/build
#            SVN_REPOSITORY https://niftyseg.svn.sourceforge.net/svnroot/niftyseg/
#            UPDATE_COMMAND ""
#            CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/niftyseg
#            BUILD_COMMAND ${CMAKE_BUILD_TOOL}
#            INSTALL_COMMAND ${CMAKE_BUILD_TOOL} install
#            )
#SUBDIRS(nifty_seg_bindings)
#ENDIF(NiftyRec_INSTALL_NIFTYSEG)
SUBDIRS(seg-lib)

#-----------------------------------------------------------------------------
IF(NiftyRec_DO_PACKAGE)
        IF(UNIX)
            IF(APPLE)
                SET(CPACK_GENERATOR "TGZ;PackageMaker")
            ELSE(APPLE)
                SET(CPACK_GENERATOR "TGZ;DEB")
            ENDIF(APPLE)
        ENDIF(UNIX)

        SET(CPACK_PACKAGE_NAME "NiftyRec")
        SET(CPACK_PACKAGE_VENDOR "CMIC-UCL")
        SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NiftyRec - Tomographic Reconstruction")
        SET(CPACK_PACKAGE_VERSION_MAJOR ${Nifty_Rec_VERSION_MAJOR})
        SET(CPACK_PACKAGE_VERSION_MINOR ${Nifty_Rec_VERSION_MINOR})
        SET(CPACK_PACKAGE_INSTALL_DIRECTORY "NiftyRec")
        SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
        SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
        #SET(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/images/niftyrec_logo.ico")
        #SET(CPACK_NSIS_MUI_ICON "${CMAKE_CURRENT_SOURCE_DIR}/images/niftyrec_logo.ico")
        #SET(CPACK_NSIS_MUI_UNIICON "${CMAKE_CURRENT_SOURCE_DIR}/images/niftyrec_logo.ico")

        IF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
          SET (CMAKE_SYSTEM_PROCESSOR "x86")
        ENDIF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
        IF(NOT DEFINED CPACK_SYSTEM_NAME)
          SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
        ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)
        IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
          IF(CMAKE_CL_64)
            SET(CPACK_SYSTEM_NAME Win64-${CMAKE_SYSTEM_PROCESSOR})
          ELSE(CMAKE_CL_64)
            SET(CPACK_SYSTEM_NAME Win32-${CMAKE_SYSTEM_PROCESSOR})
          ENDIF(CMAKE_CL_64)
        ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)

        IF(${CPACK_SYSTEM_NAME} MATCHES Darwin AND CMAKE_OSX_ARCHITECTURES)
          list(LENGTH CMAKE_OSX_ARCHITECTURES _length)
          IF(_length GREATER 1)
            SET(CPACK_SYSTEM_NAME Darwin-Universal)
          ELSE(_length GREATER 1)
            SET(CPACK_SYSTEM_NAME Darwin-${CMAKE_OSX_ARCHITECTURES})
          ENDIF(_length GREATER 1)
        ENDIF(${CPACK_SYSTEM_NAME} MATCHES Darwin AND CMAKE_OSX_ARCHITECTURES)

        #Deb
        SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Stefano Pedemonte s.pedemonte@cs.ucl.ac.uk")
        #SET(CPACK_DEBIAN_PACKAGE_DEPENDS "niftyreg (>= 1.2.1)")

        #Nsis
        SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\NiftyRec.exe")
        SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} NiftyRec")
        SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\cmic.cs.ucl.ac.uk")
        SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\cmic.cs.ucl.ac.uk")
        SET(CPACK_NSIS_CONTACT "s.pedemonte@cs.ucl.ac.uk")
        SET(CPACK_NSIS_MODIFY_PATH ON)

        #Apple drag-and-drop
        IF(APPLE)
                #SET(CPACK_PACKAGE_INSTALL_DIRECTORY "/Library/NiftyRec/")
                SET(CPACK_BINARY_TBZ2 OFF)
                SET(CPACK_BINARY_DRAGNDROP ON)
                SET(CPACK_BINARY_PACKAGEMAKER ON)
                SET(CPACK_BINARY_STGZ OFF)
        ENDIF(APPLE)

        INCLUDE(CPack)
ENDIF(NiftyRec_DO_PACKAGE)

#-----------------------------------------------------------------------------

