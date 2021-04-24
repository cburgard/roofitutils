file(GLOB TestScripts "${CMAKE_CURRENT_SOURCE_DIR}/test/*.sh")
foreach(TestScript ${TestScripts})
  if(${ATLAS_BUILD})
    add_test(NAME ${TestName} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND bash -e ${TestScript})
  else()
    atlas_add_test(NAME ${TestName} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND bash -e ${TestScript})
  endif()
  set_tests_properties(${TestName} PROPERTIES ENVIRONMENT "PYTHONPATH=$ENV{PYTHONPATH}:${EXPORT_PYTHONPATH}")
endforeach()

file(GLOB PythonTests "${CMAKE_CURRENT_SOURCE_DIR}/test/*.py")
foreach(Test ${PythonTests})
  get_filename_component(TestName ${Test} NAME)
  if(${ATLAS_BUILD})
    add_test(NAME ${TestName} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND python ${Test})
  else()
    atlas_add_test(NAME ${TestName} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND python ${Test})
  endif()
  set_tests_properties(${TestName} PROPERTIES ENVIRONMENT "PYTHONPATH=$ENV{PYTHONPATH}:${EXPORT_PYTHONPATH}")
endforeach()

enable_testing()
