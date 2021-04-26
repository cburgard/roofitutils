#########################################################################################################
# Package: RooFitUtils ###################################################################################

cmake_minimum_required(VERSION 3.0)
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

include_directories ("${PROJECT_SOURCE_DIR}")
find_package( ROOT REQUIRED COMPONENTS Core RIO MathCore Matrix HistFactory RooFitCore RooFit Hist RooStats Minuit2 Minuit )
FOREACH(incfile ${ROOT_USE_FILE})
  include(${incfile})
ENDFOREACH()

# set some variables for easier handling
set(RooFitUtilsLinkDef Root/LinkDef.h)
file(GLOB RooFitUtilsSources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "Root/[A-Za-z]*.cxx")
file(GLOB RooFitUtilsHeaders RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "RooFitUtils/[A-Za-z]*.h")


find_package(RooFitExtensions QUIET)
if(${RooFitExtensions_FOUND})
  list(APPEND CMAKE_PREFIX_PATH $ENV{RooFitExtensions_DIR})
  include_directories(${RooFitExtensions_INCLUDE_DIRS})
  link_directories(${RooFitExtensions_LIBRARY_DIRS}) 
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DHAS_ROOSTARMOMENTMORPH")
else()
  message("disabling RooStarMomentmorph integration")
endif()

if(ATLAS_BUILD OR NOT STANDALONE_BUILD)
  include(cmake/ASG.cmake)
endif()

if(STANDALONE_BUILD OR NOT ATLAS_BUILD)
  include(cmake/standalone.cmake)
endif()

include(cmake/testing.cmake)

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )
file(APPEND ${SETUP} "ulimit -S -s unlimited\n" )
file(APPEND ${SETUP} "export PATH=\${PATH}:${CMAKE_CURRENT_SOURCE_DIR}/scripts:${CMAKE_CURRENT_SOURCE_DIR}/share\n")
file(APPEND ${SETUP} "export PYTHONPATH=\${PYTHONPATH}:${EXPORT_PYTHONPATH}\n")
file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${EXPORT_LD_LIBRARY_PATH}\n")
file(APPEND ${SETUP} "export DYLD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${EXPORT_LD_LIBRARY_PATH}\n")
file(APPEND ${SETUP} "export ROOT_INCLUDE_PATH=\${ROOT_INCLUDE_PATH}:${EXPORT_ROOT_INCUDE_PATH}\n")

# general post-compile installation of the package
# required for other packages to discover this one via CMake
install(
  TARGETS RooFitUtils
  EXPORT RooFitUtilsConfig
  DESTINATION lib
  )
