#########################################################################################################
# Package: RooFitUtils ###################################################################################

cmake_minimum_required(VERSION 3.5)
project( RooFitUtils )

# figure out if we're the top-level package
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(HAS_PARENT 0)
else()
  set(HAS_PARENT 1)
endif()

set(CMAKE_BUILD_TYPE RelWithDebInfo)

IF(${HAS_PARENT})
  message("announcing RooFitUtils")
  set(HAS_RooFitUtils 1 PARENT_SCOPE)
ENDIF()

# add the current source directory to the include path
include_directories ("${PROJECT_SOURCE_DIR}")

# find ROOT (which we depend on) including all required components
find_package( ROOT REQUIRED COMPONENTS Core RIO MathCore Matrix HistFactory RooFitCore RooFit Hist RooStats Minuit2 Minuit OPTIONAL_COMPONENTS RooFitHS3 )
# include the ROOT_USE_FILE to get access to some of the ROOT-specific CMake magic
FOREACH(incfile ${ROOT_USE_FILE})
  include(${incfile})
ENDFOREACH()

# set some variables for easier handling
set(RooFitUtilsLinkDef Root/LinkDef.h)
# we glob for our source and header files so we don't have to update the lists by hand in case anything ever changes
file(GLOB RooFitUtilsSources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "Root/[A-Za-z]*.cxx")
file(GLOB RooFitUtilsHeaders RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "RooFitUtils/[A-Za-z]*.h")

# we check if RooFitExtensions is there
# if it is, we can enable some extra functionality
find_package(RooFitExtensions QUIET NO_DEFAULT_PATH) 
if(${RooFitExtensions_FOUND})
  message("found RooFitExtensions")
  list(APPEND CMAKE_PREFIX_PATH $ENV{RooFitExtensions_DIR})
  include_directories(${RooFitExtensions_INCLUDE_DIRS})
  link_directories(${RooFitExtensions_LIBRARY_DIRS})
  message(${RooFitExtensions_LIBRARY_DIRS}) 
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DHAS_ROOSTARMOMENTMORPH")
else()
  message("disabling RooStarMomentmorph integration")
endif()

# we check if we are inside an ATLAS environment
# if yes, we include the ASG CMake file
if(ATLAS_BUILD OR NOT STANDALONE_BUILD)
  include(cmake/ASG.cmake)
endif()

# if we build standalone (=without ATLAS), we use the standalone CMake file
if(STANDALONE_BUILD OR NOT ATLAS_BUILD)
  include(cmake/standalone.cmake)
endif()

# add the extra CMake snippet to register the test cases
include(cmake/testing.cmake)

# create the setup.sh script to source all the paths required by this package
set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )
# this is a neat little hack that will allow to get around some memory allocation issues with RooWorkspaces
# not really related to the package, but nice to have anyway
file(APPEND ${SETUP} "ulimit -S -s unlimited\n" )
# add PATH to be able to execute our executables
file(APPEND ${SETUP} "export PATH=\${PATH}:${CMAKE_CURRENT_SOURCE_DIR}/scripts:${CMAKE_CURRENT_SOURCE_DIR}/share\n")
# add PYTHONPATH to allow to import our python modules
file(APPEND ${SETUP} "export PYTHONPATH=\${PYTHONPATH}:${EXPORT_PYTHONPATH}\n")
# add LD_LIBRARY_PATH to allow linking to our .so library files
file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${EXPORT_LD_LIBRARY_PATH}\n")
# same thing but for MacOS
file(APPEND ${SETUP} "export DYLD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${EXPORT_LD_LIBRARY_PATH}\n")
# add some extra stuff to the ROOT include path, which is sometimes used by ROOT dictionary generation
file(APPEND ${SETUP} "export ROOT_INCLUDE_PATH=\${ROOT_INCLUDE_PATH}:${EXPORT_ROOT_INCUDE_PATH}\n")

# general post-compile installation of the package
# required for other packages to discover this one via CMake
install(
  TARGETS RooFitUtils
  EXPORT RooFitUtilsConfig
  DESTINATION lib
  )
