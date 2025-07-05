# Source directory CMakeLists.txt for TagDust2

# Include the configured header
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Define that we have a config header
add_compile_definitions(HAVE_CONFIG_H=1)

# Define common source groups
set(COMMON_SOURCES
    kslib.c
    nuc_code.c 
    misc.c
    interface.c
)

set(IO_SOURCES
    io.c
)

# Main tagdust executable
add_executable(tagdust
    ${COMMON_SOURCES}
    ${IO_SOURCES}
    main.c
    barcode_hmm.c
    calibrateQ.c
    test_architectures.c
)

target_compile_definitions(tagdust PRIVATE TAGDUST)
target_link_libraries(tagdust ${CMAKE_THREAD_LIBS_INIT} ${MATH_LIBRARY})

# simreads executable
add_executable(simreads
    ${COMMON_SOURCES}
    simulate_reads.c
)

target_compile_definitions(simreads PRIVATE SIMREADS)
target_link_libraries(simreads ${MATH_LIBRARY})

# evalres executable  
add_executable(evalres
    ${COMMON_SOURCES}
    ${IO_SOURCES}
    evalresults.c
)

target_compile_definitions(evalres PRIVATE EVALRES)
target_link_libraries(evalres ${MATH_LIBRARY})

# merge executable
add_executable(merge
    ${COMMON_SOURCES}
    ${IO_SOURCES}
    merge.c
)

target_compile_definitions(merge PRIVATE MERGE)
target_link_libraries(merge ${CMAKE_THREAD_LIBS_INIT} ${MATH_LIBRARY})

# rename_qiime executable
add_executable(rename_qiime
    rename_illumina_demultiplexed_for_qiime.c
    nuc_code.c
    kslib.c
    interface.c
    misc.c
    io.c
)

target_compile_definitions(rename_qiime PRIVATE RENAME)
target_link_libraries(rename_qiime ${MATH_LIBRARY})

# Test executables (only built when testing is enabled)
if(BUILD_TESTS)
    # Main test program
    add_executable(tagdustiotest
        io.c
        misc.c
        interface.c
        kslib.c
    )
    
    target_compile_definitions(tagdustiotest PRIVATE UTEST TAGDUST)
    target_link_libraries(tagdustiotest ${CMAKE_THREAD_LIBS_INIT} ${MATH_LIBRARY})

    # Reproducible test versions with different random seed
    add_executable(tagdust_rtest
        ${COMMON_SOURCES}
        ${IO_SOURCES}
        main.c
        barcode_hmm.c
        calibrateQ.c
        test_architectures.c
    )
    
    target_compile_definitions(tagdust_rtest PRIVATE RTEST TAGDUST)
    target_link_libraries(tagdust_rtest ${CMAKE_THREAD_LIBS_INIT} ${MATH_LIBRARY})

    add_executable(simreads_rtest
        ${COMMON_SOURCES}
        simulate_reads.c
    )
    
    target_compile_definitions(simreads_rtest PRIVATE RTEST SIMREADS)
    target_link_libraries(simreads_rtest ${MATH_LIBRARY})

    add_executable(evalres_rtest
        ${COMMON_SOURCES}
        ${IO_SOURCES}
        evalresults.c
    )
    
    target_compile_definitions(evalres_rtest PRIVATE RTEST EVALRES)
    target_link_libraries(evalres_rtest ${MATH_LIBRARY})

    # Add basic unit test
    add_test(NAME tagdustiotest COMMAND tagdustiotest)
    
    # Set test environment variables
    set_tests_properties(tagdustiotest PROPERTIES
        ENVIRONMENT "devdir=${CMAKE_SOURCE_DIR}/dev"
    )
    
    if(ENABLE_VALGRIND_TESTS AND VALGRIND_EXECUTABLE)
        set_tests_properties(tagdustiotest PROPERTIES
            COMMAND "${VALGRIND_COMMAND} $<TARGET_FILE:tagdustiotest>"
        )
    endif()
endif()

# Installation of binaries
install(TARGETS 
    tagdust 
    simreads 
    evalres 
    merge 
    rename_qiime
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install headers for development (optional)
install(FILES
    barcode_hmm.h
    interface.h
    io.h
    kslib.h
    misc.h
    nuc_code.h
    tagdust2.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tagdust2
    OPTIONAL
)