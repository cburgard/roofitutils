file(GLOB TestScripts "${CMAKE_CURRENT_SOURCE_DIR}/test/*.sh")

foreach(TestScript ${TestScripts})
  get_filename_component(TestName ${TestScript} NAME)
  add_test(NAME ${TestName} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND bash -e ${TestScript})
  set_tests_properties(${TestName} PROPERTIES ENVIRONMENT "PYTHONPATH=$ENV{PYTHONPATH}:${EXPORT_PYTHONPATH};LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${EXPORT_LD_LIBRARY_PATH}")
endforeach()

file(GLOB PythonTests "${CMAKE_CURRENT_SOURCE_DIR}/test/*.py")
foreach(Test ${PythonTests})
  get_filename_component(TestName ${Test} NAME)
  add_test(NAME ${TestName} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMAND python ${Test})
  set_tests_properties(${TestName} PROPERTIES ENVIRONMENT "PYTHONPATH=$ENV{PYTHONPATH}:${EXPORT_PYTHONPATH};LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${EXPORT_LD_LIBRARY_PATH}")
endforeach()

enable_testing()
