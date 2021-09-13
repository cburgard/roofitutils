# register all the files and directories
include_directories ("${PROJECT_SOURCE_DIR}" "${ROOT_INCLUDE_DIRS}")

# generate the dictionary source code
ROOT_GENERATE_DICTIONARY(G__RooFitUtils ${RooFitUtilsHeaders} LINKDEF ${RooFitUtilsLinkDef})

# register the shared object to include both sources and dictionaries
add_library( RooFitUtils SHARED ${RooFitUtilsSources} G__RooFitUtils.cxx)

# link everything together at the end
target_link_libraries( RooFitUtils ${ROOT_LIBRARIES} )

file(GLOB pyfiles "python/*.py")
execute_process(
  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/RooFitUtils
  )
foreach(pyfile ${pyfiles})
  execute_process(
    COMMAND ln -sf ${pyfile} ${CMAKE_CURRENT_BINARY_DIR}/RooFitUtils
    )
endforeach()

# Add all targets to the build-tree export set
export(TARGETS RooFitUtils FILE "${PROJECT_BINARY_DIR}/RooFitUtilsTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE RooFitUtils)

set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/RooFitUtils" )
set(CONF_LIBRARY_DIRS "${PROJECT_BINARY_DIR}")
set(CONF_LIBRARIES    RooFitUtils)
configure_file(RooFitUtilsConfig.cmake.in
  "${PROJECT_BINARY_DIR}/RooFitUtilsConfig.cmake" @ONLY)

# Install the RooFitUtilsConfig.cmake
install(FILES
  "${PROJECT_BINARY_DIR}/RooFitUtilsConfig.cmake"
  DESTINATION "${PROJECT_SOURCE_DIR}" COMPONENT dev
  )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libRooFitUtils_rdict.pcm DESTINATION lib)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libRooFitUtils.rootmap DESTINATION lib)


file(GLOB pyfiles "scripts/*")
foreach(pyfile ${pyfiles})
  install(FILES
    ${pyfile}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE WORLD_READ WORLD_EXECUTE GROUP_READ GROUP_EXECUTE
    DESTINATION bin
    )
endforeach()

install(FILES
  "${PROJECT_BINARY_DIR}/RooFitUtilsConfig.cmake"
  DESTINATION cmake
  )

# install the python package
install(DIRECTORY
  "${PROJECT_BINARY_DIR}/RooFitUtils"
  DESTINATION lib
  )

set(EXPORT_PYTHONPATH ${CMAKE_CURRENT_BINARY_DIR})
set(EXPORT_LD_LIBRARY_PATH ${CMAKE_CURRENT_BINARY_DIR})
set(EXPORT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
set(EXPORT_ROOT_INCLUDE_PATH ${CMAKE_BINARY_DIR}/include)

if(${RooFitExtensions_FOUND})
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads REQUIRED)
  target_link_libraries( RooFitUtils ${ROOT_LIBRARIES} ${RooFitExtensions_LIBRARIES} Threads::Threads)
endif()

