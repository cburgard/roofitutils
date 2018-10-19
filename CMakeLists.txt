#########################################################################################################
# Package: RooFitUtils ###################################################################################

cmake_minimum_required(VERSION 2.8)
project( RooFitUtils )

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(HAS_PARENT 0)
else()
  set(HAS_PARENT 1)
endif()
  

IF(${HAS_PARENT})
  message("announcing RooFitUtils")
  set(HAS_RooFitUtils 1 PARENT_SCOPE)
ENDIF()

find_package( AnalysisBase QUIET) 

IF(${AnalysisBase_FOUND})
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
ENDIF()

find_package( ROOT REQUIRED COMPONENTS Core RIO MathCore Matrix HistFactory RooFitCore RooFit Hist RooStats )
FOREACH(incfile ${ROOT_USE_FILE})
  include(${incfile})
ENDFOREACH()  

# set some variables for easier handling
set(RooFitUtilsLinkDef ${PROJECT_SOURCE_DIR}/Root/LinkDef.h)
file(GLOB RooFitUtilsSources Root/[A-Z]*.cxx)
file(GLOB RooFitUtilsHeaders RooFitUtils/[A-Z]*.h)

if(EXISTS ${TQPATH})
  # generate the definitions.h
  execute_process(
    COMMAND python ${TQPATH}/share/generateLocals.py --verbose  --output ${PROJECT_SOURCE_DIR}/Root/definitions.h --set-working-directory ${PROJECT_SOURCE_DIR} --flag-header HAS_ROOLAGRANGIANMORPHING RooLagrangianMorphing.h 
    )
else()
  execute_process(
    COMMAND touch ${PROJECT_SOURCE_DIR}/Root/definitions.h
    )
endif()

set(SETUP setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )

if(${AnalysisBase_FOUND})
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

  set(BINARY_TAG ${ATLAS_PLATFORM})

  file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/${BINARY_TAG}/lib\n")
  file(APPEND ${SETUP} "export DYLD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/${BINARY_TAG}/lib\n")
  file(APPEND ${SETUP} "export ROOT_INCLUDE_PATH=\${ROOT_INCLUDE_PATH}:${CMAKE_BINARY_DIR}/${BINARY_TAG}/include\n")

ELSE()
  # this section is a standalone compilation segment
  # it is executed when compiling outside an ASG environment

  # register all the files and directories
  include_directories ("${PROJECT_SOURCE_DIR}" "${ROOT_INCLUDE_DIRS}")

  # generate the dictionary source code
  ROOT_GENERATE_DICTIONARY(G__RooFitUtils ${RooFitUtilsHeaders} LINKDEF ${RooFitUtilsLinkDef})

  # register the shared object to include both sources and dictionaries
  add_library( RooFitUtils SHARED ${RooFitUtilsSources} G__RooFitUtils.cxx)

  # link everything together at the end
  target_link_libraries( RooFitUtils ${ROOT_LIBRARIES} )

  file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_CURRENT_BINARY_DIR}\n")
  file(APPEND ${SETUP} "export PATH=\${PATH}:${CMAKE_CURRENT_SOURCE_DIR}/scripts\n")
  
ENDIF()

# general post-compile installation of the package 
# required for other packages to discover this one via CMake
install( 
  TARGETS RooFitUtils
  EXPORT RooFitUtilsConfig
  DESTINATION cmake
  )
