# register all the files and directories
include_directories ("${PROJECT_SOURCE_DIR}" "${ROOT_INCLUDE_DIRS}")

# generate the dictionary source code
ROOT_GENERATE_DICTIONARY(G__RooFitUtils ${RooFitUtilsHeaders} LINKDEF ${RooFitUtilsLinkDef})

# register the shared object to include both sources and dictionaries
add_library( RooFitUtils SHARED ${RooFitUtilsSources} G__RooFitUtils.cxx)

# link everything together at the end
target_link_libraries( RooFitUtils ${ROOT_LIBRARIES} )

# create a folder to host all of our python submodules
execute_process(
  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/RooFitUtils
  )
# find all of our python submodules and link them to the folder above
file(GLOB pyfiles "python/*.py")
foreach(pyfile ${pyfiles})
  execute_process(
    COMMAND ln -sf ${pyfile} ${CMAKE_CURRENT_BINARY_DIR}/RooFitUtils
    )
endforeach()

# add all targets to the build-tree export set
# this is needed for other packages to find this one, in case they depend on us
export(TARGETS RooFitUtils FILE "${PROJECT_BINARY_DIR}/RooFitUtilsTargets.cmake")
export(PACKAGE RooFitUtils)

# here we create the RooFitUtils.Config.cmake file
# this is needed for other packages to find this one with FIND_PACKAGE
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/RooFitUtils" )
set(CONF_LIBRARY_DIRS "${PROJECT_BINARY_DIR}")
set(CONF_LIBRARIES    RooFitUtils)
configure_file(RooFitUtilsConfig.cmake.in
  "${PROJECT_BINARY_DIR}/RooFitUtilsConfig.cmake" @ONLY)
# here we install this file to a path where it can be found
install(FILES
  "${PROJECT_BINARY_DIR}/RooFitUtilsConfig.cmake"
  DESTINATION "${PROJECT_SOURCE_DIR}" COMPONENT dev
  )
install(FILES
  "${PROJECT_BINARY_DIR}/RooFitUtilsConfig.cmake"
  DESTINATION cmake
  )

# we also need to install the .pcm and .rootmap files to help ROOT find the classes in this package
# this is needed to get rid of some warning messages when using ROOT to read in root files when this package is loaded
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libRooFitUtils_rdict.pcm DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libRooFitUtils.rootmap DESTINATION lib)

# we also loop over all the scripts in this package and install them
file(GLOB pyfiles "scripts/*")
foreach(pyfile ${pyfiles})
  install(FILES
    ${pyfile}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE WORLD_READ WORLD_EXECUTE GROUP_READ GROUP_EXECUTE
    DESTINATION bin
    )
endforeach()

# install the python package
install(DIRECTORY
  "${PROJECT_BINARY_DIR}/RooFitUtils"
  DESTINATION lib
  )

# we set the export paths so that they can be used by the parent CMakeLists.txt
set(EXPORT_PYTHONPATH ${CMAKE_CURRENT_BINARY_DIR})
set(EXPORT_LD_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR})
set(EXPORT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
set(EXPORT_ROOT_INCLUDE_PATH ${CMAKE_BINARY_DIR}/include)

# if we have RooFitExtensions available, we need to make sure we include its dependencies
if(${RooFitExtensions_FOUND})
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  target_link_libraries( RooFitUtils ${ROOT_LIBRARIES} ${RooFitExtensions_LIBRARIES} Threads::Threads)
endif()


