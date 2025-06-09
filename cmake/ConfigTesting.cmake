set(BUILDNAME
    "${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${CMAKE_Fortran_COMPILER_ID}-${BLAS_TYPE}-${CMAKE_BUILD_TYPE}"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest own timeout
if(ENABLE_LARGE_TEST)
   set(DART_TESTING_TIMEOUT
      "18000" # 5 hours 
      CACHE STRING
      "Set timeout in seconds for every single test"
      )
else()
   set(DART_TESTING_TIMEOUT
      "1200"  # 20 minutes
      CACHE STRING
      "Set timeout in seconds for every single test"
      )
endif()

include(Tests_plasmonX)
include(CTest)
enable_testing()

