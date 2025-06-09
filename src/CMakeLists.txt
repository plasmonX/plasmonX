# Ensure correct policy settings
cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0079 NEW)

# Create a static library for shared modules
add_library(plasmonX_modules STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/bem/bem_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/bem/exact_green_function.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/wfq/wfq_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/wfq/wfqfmu_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/forcefield/fq_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/forcefield/fqfmu_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/control/control_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/solvers/algorithm_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/solvers/gmres_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/input_output/output_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/input_output/input_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/parameters/parameters_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/matrix/matrix_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/memory_manager/memory_manager_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/array_manipulation_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/float_manipulation_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tools/string_manipulation_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/rhs/rhs_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/field/field_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/target/target_module.f90
)

# Ensure module files are available
target_include_directories(plasmonX_modules PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Define source files for plasmonX and pass them to the parent scope
set(plasmonX_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/plasmonX.f90
)
set(plasmonX_ANALYSIS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/plasmonX_analysis/input_analysis_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/plasmonX_analysis/output_analysis_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/plasmonX_analysis/control_analysis_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/plasmonX_analysis/grid_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/plasmonX_analysis/plot_module.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/plasmonX_analysis/plasmonX_analysis.f90
)

# Pass variables to the parent CMakeLists
set(plasmonX_SOURCES ${plasmonX_SOURCES} PARENT_SCOPE)
set(plasmonX_ANALYSIS_SOURCES ${plasmonX_ANALYSIS_SOURCES} PARENT_SCOPE)
