cmake_minimum_required(VERSION 3.15)

project(plasmonX)
project(plasmonX_analysis)
enable_language(Fortran)


set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/math
    /usr/local/share/cmake/Modules
    )

#set(CMAKE_MODULE_PATH "/opt/cmake-3.14.0-rc1-Linux-x86_64/share/cmake-3.14/Modules/;${CMAKE_MODULE_PATH}")
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)
option(ENABLE_THREADED_MKL   "Enable OpenMP parallelization in MKL"             ON)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    OFF)
option(ENABLE_INT64          "Enable 64-bit integers"                           ON)

# Build type options
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# ============================
# Build Type and Compilation Flags
# ============================

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (Release, ReleaseWithWarnings, Debug)." FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Define compilation flags based on build type
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -funroll-loops -std=f2008")
    elseif(CMAKE_BUILD_TYPE STREQUAL "ReleaseWithWarnings")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -Wall -Wextra -Wimplicit-interface -Wconversion -std=f2008")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -Wall -Wextra -Wimplicit-interface -Wconversion -fcheck=all -finit-real=snan -finit-integer=-999999 -finit-logical=false -std=f2008")
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -xHost -funroll-loops -ipo -std=f2008")
    elseif(CMAKE_BUILD_TYPE STREQUAL "ReleaseWithWarnings")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -warn all -warn interfaces -traceback -std=f2008")
    elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -traceback -debug full -check bounds -warn all -warn interfaces -fpe0 -init=snan -init=integer=-999999 -std=f2008")
    endif()
endif()

# ========================
# Ensure OpenMP is Enabled
# ========================

if(ENABLE_OMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    if(OpenMP_Fortran_FOUND)
        message(STATUS "OpenMP Fortran found: ${OpenMP_Fortran_FLAGS}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")
    else()
        message(FATAL_ERROR "OpenMP Fortran NOT FOUND! Check compiler or flags.")
    endif()
else()
    message(WARNING "OpenMP is DISABLED! Re-run cmake with -DENABLE_OMP=ON if needed.")
endif()


include(ConfigTesting)

# Math detection
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)

if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#
# Compile.
#
# Create the executables using the source files from src/
add_subdirectory(src) 

# Create the executables
add_executable(plasmonX ${plasmonX_SOURCES})
add_executable(plasmonX_analysis ${plasmonX_ANALYSIS_SOURCES})

# Ensure correct module directories
target_include_directories(plasmonX PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_include_directories(plasmonX_analysis PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)

# Link the shared module library
target_link_libraries(plasmonX PRIVATE plasmonX_modules)
target_link_libraries(plasmonX_analysis PRIVATE plasmonX_modules)

# ========================
# Link OpenMP and BLAS/LAPACK
# ========================
if(ENABLE_OMP)
    foreach(EXE plasmonX plasmonX_analysis)
        target_link_libraries(${EXE} 
            PUBLIC OpenMP::OpenMP_Fortran 
            PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
        )
    endforeach()
else()
    foreach(EXE plasmonX plasmonX_analysis)
        target_link_libraries(${EXE} PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    endforeach()
endif()


find_package(Python3 REQUIRED COMPONENTS Interpreter) 
  
file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})

#runtest
if(EXISTS ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest)
   file(COPY ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest DESTINATION ${CMAKE_BINARY_DIR}/tests)
endif()

# Assicura che il submodule GEOM esista
if(NOT EXISTS "${CMAKE_SOURCE_DIR}/submodules/geom")
    message(FATAL_ERROR "GEOM submodule not found! Run: git submodule update --init")
endif()

# Assicura che Python trovi GEOM nel PYTHONPATH
set(ENV{PYTHONPATH} "${CMAKE_SOURCE_DIR}/submodules/geom:$ENV{PYTHONPATH}")

# Percorso per i requirements e i test
set(PYTHON_REQUIREMENTS "${CMAKE_SOURCE_DIR}/python/requirements.txt")
set(GEOM_TESTS_PATH "${CMAKE_SOURCE_DIR}/submodules/geom/geom/tests")

# Controlla se tutti i requirements sono già installati
if(EXISTS "${PYTHON_REQUIREMENTS}")
    message(STATUS "Checking GEOM dependencies...")

    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip check
        RESULT_VARIABLE REQUIREMENTS_OK
        OUTPUT_VARIABLE REQUIREMENTS_OUTPUT
        ERROR_VARIABLE REQUIREMENTS_OUTPUT
    )

    if (REQUIREMENTS_OK)
        message(STATUS "Installing missing GEOM dependencies globally...")

        execute_process(
            COMMAND ${Python3_EXECUTABLE} -m pip install -q --disable-pip-version-check --no-cache-dir -r ${PYTHON_REQUIREMENTS}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/submodules/geom
            RESULT_VARIABLE REQUIREMENTS_FAILED
        )

        if (REQUIREMENTS_FAILED)
            message(FATAL_ERROR "Failed to install GEOM dependencies!")
        endif()
    else()
        message(STATUS "GEOM dependencies are already installed.")
    endif()
else()
    message(WARNING "requirements.txt not found in ${GEOM_REQUIREMENTS}. Skipping dependency installation.")
endif()

# Assicuriamoci che i test di GEOM esistano prima di procedere
if(NOT EXISTS "${GEOM_TESTS_PATH}/run_all_tests.sh")
    message(FATAL_ERROR "GEOM tests not found in ${GEOM_TESTS_PATH}.")
endif()

# File che memorizza il commit e lo stato dei test di GEOM
set(GEOM_TEST_STATUS_FILE "${CMAKE_BINARY_DIR}/CMake_GEOM_Tests.txt")

# Ottieni l'ultimo commit del submodule GEOM
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/submodules/geom
    OUTPUT_VARIABLE GEOM_CURRENT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Variabile per decidere se eseguire i test
set(SHOULD_RUN_GEOM_TESTS ON)

# Se il file di stato esiste, controlliamo il commit e lo stato dei test
if(EXISTS "${GEOM_TEST_STATUS_FILE}")
    file(STRINGS "${GEOM_TEST_STATUS_FILE}" GEOM_TEST_STATUS_LIST)
    list(GET GEOM_TEST_STATUS_LIST 0 SAVED_GEOM_COMMIT)
    list(GET GEOM_TEST_STATUS_LIST 1 SAVED_TEST_RESULT)

    # Rimuoviamo spazi e newline extra dalle variabili
    string(STRIP "${SAVED_GEOM_COMMIT}" SAVED_GEOM_COMMIT)
    string(STRIP "${GEOM_CURRENT_COMMIT}" GEOM_CURRENT_COMMIT)
    string(STRIP "${SAVED_TEST_RESULT}" SAVED_TEST_RESULT)
    
    # Controlliamo se il commit è uguale e se la stringa di test è corretta
    if("${SAVED_GEOM_COMMIT}" STREQUAL "${GEOM_CURRENT_COMMIT}" AND "${SAVED_TEST_RESULT}" MATCHES "TESTS_PASSED")
        message(STATUS "GEOM tests already run and passed for commit ${GEOM_CURRENT_COMMIT}. Skipping tests.")
        set(SHOULD_RUN_GEOM_TESTS OFF)
    endif()
endif()

# Se dobbiamo eseguire i test, li lanciamo
if(SHOULD_RUN_GEOM_TESTS)
    message(STATUS "Running GEOM tests...")

    execute_process(
        COMMAND bash -c "PYTHONPATH=${CMAKE_SOURCE_DIR}/submodules/geom:\$PYTHONPATH ${GEOM_TESTS_PATH}/run_all_tests.sh > ${CMAKE_BINARY_DIR}/geom_test_log.txt 2>&1"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        RESULT_VARIABLE TEST_RESULT
    )

    # Controlliamo se i test hanno avuto errori
    file(STRINGS "${CMAKE_BINARY_DIR}/geom_test_log.txt" TEST_LOG_CONTENT)
    string(FIND "${TEST_LOG_CONTENT}" "ERROR" ERROR_FOUND)

    if (ERROR_FOUND GREATER -1)
        message(FATAL_ERROR "GEOM tests failed! Check ${CMAKE_BINARY_DIR}/geom_test_log.txt")
    else()
        message(STATUS "GEOM tests passed!")
        # Salva il commit corrente e il risultato positivo dei test nel file di stato
        file(WRITE "${GEOM_TEST_STATUS_FILE}" "${GEOM_CURRENT_COMMIT}\n")
        file(APPEND "${GEOM_TEST_STATUS_FILE}" "TESTS_PASSED\n")
    endif()
endif()



#get git version
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT GIT_BRANCH)
    set(GIT_BRANCH "unknown-branch")
endif()


# Definition for pointers in plasmonX.py
set(BUILD_PATH ${CMAKE_BINARY_DIR})
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})
#time of compilation
execute_process(COMMAND date +%Z OUTPUT_VARIABLE SYSTEM_TIMEZONE OUTPUT_STRIP_TRAILING_WHITESPACE)
string(TIMESTAMP CONFIGURATION_DATE "%Y-%m-%d %H:%M:%S ${SYSTEM_TIMEZONE}")
set(CONFIGURATION_DATE "${CONFIGURATION_DATE}")
set(GIT_BRANCH ${GIT_BRANCH})  
set(FORTRAN_COMPILER "${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
#LAPACK and BLAS
function(detect_library_type library_list result_var)
    if("${library_list}" MATCHES "mkl")
        set(${result_var} "MKL" PARENT_SCOPE)
    elseif("${library_list}" MATCHES "openblas")
        set(${result_var} "OPENBLAS" PARENT_SCOPE)
    elseif("${library_list}" MATCHES "atlas")
        set(${result_var} "ATLAS" PARENT_SCOPE)
    else()
        set(${result_var} "NATIVE" PARENT_SCOPE)
    endif()
endfunction()
if(LAPACK_FOUND AND BLAS_FOUND)
    detect_library_type("${LAPACK_LIBRARIES}" LAPACK_TYPE)
    detect_library_type("${BLAS_LIBRARIES}" BLAS_TYPE)
else()
    set(LAPACK_TYPE "unknown")
    set(BLAS_TYPE "unknown")
endif()
set(INT64_STATUS ${ENABLE_INT64})
set(OMP_STATUS ${ENABLE_OMP})

# Copy and configure plasmonX.py file
configure_file(
    ${CMAKE_SOURCE_DIR}/python/input_/plasmonX.py # File sorgente
    ${CMAKE_BINARY_DIR}/plasmonX.py               # File destinazione
    @ONLY                                       # Sostituisci solo i segnaposto
)

# Copy and configure plasmonX_analysis.py file
configure_file(
    ${CMAKE_SOURCE_DIR}/python/analysis/plasmonX_analysis.py # File sorgente
    ${CMAKE_BINARY_DIR}/plasmonX_analysis.py    # File destinazione
    @ONLY                                     # Sostituisci solo i segnaposto
)
