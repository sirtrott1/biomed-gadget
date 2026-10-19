add_library(gadget_source_formats STATIC bindingdb.cpp pubchem.cpp gtex.cpp)
target_link_libraries(gadget_source_formats PUBLIC gadget::core)
target_include_directories(gadget_source_formats PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(gadget_source_formats PRIVATE -Wall -Wextra -Wpedantic -Werror)

if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  foreach(reader bindingdb_bulk pubchem_assay gtex_source)
    add_executable(gadget_${reader}_tests "${PROJECT_SOURCE_DIR}/tests/${reader}.cpp")
    target_link_libraries(gadget_${reader}_tests PRIVATE gadget_source_formats)
    target_compile_options(gadget_${reader}_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)
    add_test(NAME gadget_${reader}_contracts COMMAND gadget_${reader}_tests)
    set_tests_properties(gadget_${reader}_contracts PROPERTIES TIMEOUT 60)
  endforeach()
  add_test(NAME gadget_source_formats_reference_contracts
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/source_formats.py"
      "$<TARGET_FILE:gadget_bindingdb_bulk_tests>" "$<TARGET_FILE:gadget_pubchem_assay_tests>"
      "$<TARGET_FILE:gadget_gtex_source_tests>")
  set_tests_properties(gadget_source_formats_reference_contracts PROPERTIES TIMEOUT 60)
endif()
