find_package( AnalysisBase QUIET)

IF(${AnalysisBase_FOUND})
  set(STANDALONE_BUILD 0)
  set(ATLAS_BUILD 1)  
  
  # Set up the usage of CTest:
  IF(NOT ${HAS_PARENT})
    atlas_ctest_setup() # Set up the project:
    atlas_project( RooFitUtils 1.0.0
      USE AnalysisBase ${AnalysisBase_VERSION} )

    # Generate an environment setup script:
    lcg_generate_env( SH_FILE ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh )
    install( FILES ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh DESTINATION . )

    # Set up the usage of CPack:
    set(CMAKE_INSTALL_PREFIX /InstallArea/x86_64-slc6-gcc49-opt)
    atlas_cpack_setup()
  ENDIF()

  # register this as a package to ASG
  atlas_subdir( RooFitUtils )

  atlas_platform_id( BINARY_TAG )

  # this section reflects the standard ASG way of configuring CMake
  # it is executed when compiling within an ASG environment
  find_package( GTest )
  set(CMAKE_INSTALL_PREFIX /InstallArea/x86_64-slc6-gcc62-opt)
  atlas_add_root_dictionary( RooFitUtils RooFitUtilsCintDict
    ROOT_HEADERS ${RooFitUtilsHeaders} ${RooFitUtilsLinkDef}
    EXTERNAL_PACKAGES ROOT )
  atlas_add_library( RooFitUtils
    ${RooFitUtilsHeaders} ${RooFitUtilsSources} ${RooFitUtilsCintDict}
    PUBLIC_HEADERS RooFitUtils
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES}
    )

  find_package(PythonInterp REQUIRED)
  atlas_install_python_modules( python/* )

  set(EXPORT_PYTHONPATH ${CMAKE_BINARY_DIR}/${BINARY_TAG}/python)
  set(EXPORT_LD_LIBRARY_PATH ${CMAKE_BINARY_DIR}/${BINARY_TAG}/lib)
  set(EXPORT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
  set(EXPORT_ROOT_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${BINARY_TAG}/include)
else()
  if(ATLAS_BUILD)
    message(FATAL_ERROR "ATLAS_BUILD set, but cannot find ASG release")
  else()
    set(ATLAS_BUILD 0)
  endif()
endif()
