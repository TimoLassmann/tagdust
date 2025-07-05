# Test directory CMakeLists.txt for TagDust2

# Find bash for running shell script tests
find_program(BASH_PROGRAM bash)

if(NOT BASH_PROGRAM)
    message(WARNING "Bash not found - shell script tests will be skipped")
    return()
endif()

# Set test environment variables
set(TEST_ENV_VARS
    "devdir=${CMAKE_CURRENT_SOURCE_DIR}"
    "PATH=${CMAKE_BINARY_DIR}/src:$ENV{PATH}"
)

if(ENABLE_VALGRIND_TESTS AND VALGRIND_EXECUTABLE)
    list(APPEND TEST_ENV_VARS "valparam=${VALGRIND_COMMAND}")
else()
    list(APPEND TEST_ENV_VARS "valparam=")
endif()

# Add shell script tests
set(SHELL_TESTS
    sanity_test.sh
    bar_read_test.sh
    casava_test.sh
)

# Create modified test scripts that use the correct executable paths
foreach(test_script ${SHELL_TESTS})
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_script}")
        get_filename_component(test_name ${test_script} NAME_WE)
        set(modified_script "${CMAKE_CURRENT_BINARY_DIR}/${test_script}")
        
        # Read the original script and modify the paths
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/${test_script}" script_content)
        
        # Replace ../src/ with the correct build path
        string(REPLACE "../src/" "${CMAKE_BINARY_DIR}/src/" script_content "${script_content}")
        
        # Write the modified script
        file(WRITE "${modified_script}" "${script_content}")
        
        # Make it executable
        file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${test_script}" 
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
             FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
        
        # Overwrite with our modified version
        file(WRITE "${modified_script}" "${script_content}")
        
        add_test(
            NAME ${test_name}
            COMMAND ${BASH_PROGRAM} ${modified_script}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        
        set_tests_properties(${test_name} PROPERTIES
            ENVIRONMENT "${TEST_ENV_VARS}"
            TIMEOUT 300  # 5 minutes timeout per test
        )
        
        # Mark as an integration test
        set_tests_properties(${test_name} PROPERTIES
            LABELS "integration"
        )
    else()
        message(STATUS "Test script ${test_script} not found - skipping")
    endif()
endforeach()

# Handle other test scripts that might exist
set(OTHER_TESTS
    edge_case_test.sh
    cli_test.sh
    bmp_test.sh
    arch_test.sh
    format_test.sh
    thread_test.sh
)

foreach(test_script ${OTHER_TESTS})
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_script}")
        get_filename_component(test_name ${test_script} NAME_WE)
        set(modified_script "${CMAKE_CURRENT_BINARY_DIR}/${test_script}")
        
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/${test_script}" script_content)
        string(REPLACE "../src/" "${CMAKE_BINARY_DIR}/src/" script_content "${script_content}")
        file(WRITE "${modified_script}" "${script_content}")
        
        file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${test_script}" 
             DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
             FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
        file(WRITE "${modified_script}" "${script_content}")
        
        add_test(
            NAME ${test_name}
            COMMAND ${BASH_PROGRAM} ${modified_script}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        
        set_tests_properties(${test_name} PROPERTIES
            ENVIRONMENT "${TEST_ENV_VARS}"
            TIMEOUT 300
            LABELS "integration"
        )
    else()
        message(STATUS "Test script ${test_script} not found - skipping")
    endif()
endforeach()

# Copy test data and architecture files for testing
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/README.txt"
    "${CMAKE_CURRENT_BINARY_DIR}/README.txt"
    COPYONLY
)

# Install test files and data (optional)
install(FILES
    README.txt
    EDITTAG_3nt_ed_1.txt
    EDITTAG_3nt_ed_2.txt
    EDITTAG_4nt_ed_2.txt
    EDITTAG_4nt_ed_3.txt
    EDITTAG_6nt_ed_3.txt
    EDITTAG_6nt_ed_4.txt
    casava_arch.txt
    barread1_tagdust_results_gold.txt
    barread2_tagdust_results_gold.txt
    barread_paired_tagdust_results_gold.txt
    read_paired_tagdust_results_gold.txt
    casava_out_BC_TTAGGC_READ1_gold.txt
    casava_out_BC_TTAGGC_READ2_gold.txt
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/tagdust2/test_data
    OPTIONAL
)

install(FILES
    casava_read1.fastq.gz
    casava_read2.fastq.gz
    casava_read3.fastq.gz
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/tagdust2/test_data
    OPTIONAL
)

install(PROGRAMS
    ${SHELL_TESTS}
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/tagdust2/tests
    OPTIONAL
)