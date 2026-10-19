add_library(gadget_source_formats STATIC bindingdb.cpp pubchem.cpp gtex.cpp pubchem_model.cpp pubchem_evidence.cpp)
target_link_libraries(gadget_source_formats PUBLIC gadget::core)
target_include_directories(gadget_source_formats PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(gadget_source_formats PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(gadget_source source.cpp)
target_link_libraries(gadget_source PRIVATE gadget_source_formats)
target_compile_options(gadget_source PRIVATE -Wall -Wextra -Wpedantic -Werror)

if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  set(source_test_data_dir "${PROJECT_BINARY_DIR}/test-data")
  foreach(native pubchem_native_evidence feature_context)
    add_executable(gadget_${native}_tests "${PROJECT_SOURCE_DIR}/tests/${native}.cpp")
    target_link_libraries(gadget_${native}_tests PRIVATE gadget_source_formats)
    target_compile_options(gadget_${native}_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)
    add_test(NAME gadget_${native}_contracts COMMAND gadget_${native}_tests "${source_test_data_dir}/${native}")
    set_tests_properties(gadget_${native}_contracts PROPERTIES TIMEOUT 120)
  endforeach()
  set_tests_properties(gadget_feature_context_contracts PROPERTIES FIXTURES_SETUP gadget_typed_evidence_inputs)
  add_test(NAME gadget_pubchem_source_cli_contracts COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/pubchem_source_cli.py"
      "$<TARGET_FILE:gadget_source>" "${source_test_data_dir}/feature_context"
      "$<TARGET_FILE:gadget>" "$<TARGET_FILE:gadget_result_binding_tests>")
  set_tests_properties(gadget_pubchem_source_cli_contracts PROPERTIES TIMEOUT 120 FIXTURES_REQUIRED gadget_typed_evidence_inputs)
  option(GADGET_TEST_NATIVE_PYTHON "Test source intake with installed native Python providers" OFF)
  if(GADGET_TEST_NATIVE_PYTHON)
    execute_process(COMMAND "${Python3_EXECUTABLE}" -I -c "import rexgraph, rexgraph_cpp, rcdb, rcql, sqlalchemy"
      RESULT_VARIABLE gadget_native_python_status ERROR_VARIABLE gadget_native_python_error)
    if(NOT gadget_native_python_status EQUAL 0)
      message(FATAL_ERROR "Native Python source tests require Core, bridge, RCDB, RCQL and SQLAlchemy: ${gadget_native_python_error}")
    endif()
    add_test(NAME gadget_source_native_contracts COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/source_native.py"
        "$<TARGET_FILE:gadget_source>" "$<TARGET_FILE:gadget>" "$<TARGET_FILE:gadget_native_field_import>" "$<TARGET_FILE:gadget_result_binding_tests>"
        "${source_test_data_dir}/source-native")
    set_tests_properties(gadget_source_native_contracts PROPERTIES TIMEOUT 180 LABELS native-python)
  endif()

  foreach(reader bindingdb_bulk pubchem_assay gtex_source)
    add_executable(gadget_${reader}_tests "${PROJECT_SOURCE_DIR}/tests/${reader}.cpp")
    target_link_libraries(gadget_${reader}_tests PRIVATE gadget_source_formats)
    target_compile_options(gadget_${reader}_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)
    add_test(NAME gadget_${reader}_contracts COMMAND gadget_${reader}_tests)
    set_tests_properties(gadget_${reader}_contracts PROPERTIES TIMEOUT 60)
  endforeach()
  add_test(NAME gadget_bindingdb_bulk_reference_contracts COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/bindingdb_bulk_reference_contracts.py" "$<TARGET_FILE:gadget_bindingdb_bulk_tests>")
  set_tests_properties(gadget_bindingdb_bulk_reference_contracts PROPERTIES TIMEOUT 60)
  add_test(NAME gadget_source_formats_reference_contracts
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/source_formats.py"
      "$<TARGET_FILE:gadget_bindingdb_bulk_tests>" "$<TARGET_FILE:gadget_pubchem_assay_tests>"
      "$<TARGET_FILE:gadget_gtex_source_tests>")
  set_tests_properties(gadget_source_formats_reference_contracts PROPERTIES TIMEOUT 60)
endif()
