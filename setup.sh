#!/bin/bash

# Colors for output
GREEN='\033[1;32m'
RED='\033[1;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No color

# Default values
build_dir="build"
build_type="Release"
enable_omp=ON
fortran_compiler=""

# Help function
usage() {
    echo "Usage: $0 [-b build_dir] [-t build_type] [-fc fortran_compiler] [-omp yes|no]"
    echo "  -b    Build directory (default: build)"
    echo "  -bt   Build type: Release, Debug, Warning (default: Release)"
    echo "  -fc   Fortran compiler (supported: gfortran)"
    echo "  -omp  Enable OpenMP: yes (default) or no"
    exit 1
}

# Parse input arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--build)
            build_dir="$2"
            shift 2
            ;;
        -bt|--build_type)
            build_type="$2"
            shift 2
            ;;
        -fc|--fortrancompiler)
            fortran_compiler="$2"
            shift 2
            ;;
        -omp|--omp)
            omp_value=$(echo "$2" | tr '[:upper:]' '[:lower:]')
            if [[ "$omp_value" == "yes" ]]; then
                enable_omp=ON
            elif [[ "$omp_value" == "no" ]]; then
                enable_omp=OFF
            else
                echo -e "${RED}Invalid value for -omp: use 'yes' or 'no'.${NC}"
                usage
            fi
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            usage
            ;;
    esac
done

# Normalize build type
build_type_lower=$(echo "$build_type" | tr '[:upper:]' '[:lower:]')

case "$build_type_lower" in
    release)
        build_type="Release"
        ;;
    debug)
        build_type="Debug"
        ;;
    warning)
        build_type="ReleaseWithWarnings"
        ;;
    *)
        echo -e "${RED}Invalid build type: choose 'Release', 'Debug', or 'Warning'.${NC}"
        usage
        ;;
esac

# Set Fortran compiler
if [[ -n "$fortran_compiler" ]]; then
    if ! command -v "$fortran_compiler" &> /dev/null; then
        echo -e "${RED}Specified Fortran compiler '$fortran_compiler' not found in PATH.${NC}"
        exit 1
    fi
    export FC="$fortran_compiler"
else
    if command -v gfortran &> /dev/null; then
        export FC=$(which gfortran)
    else
        echo -e "${YELLOW}No Fortran compiler specified and 'gfortran' not found.${NC}"
        echo -e "${YELLOW}Please note: the code is tested with 'gfortran'. Use other compilers with caution.${NC}"
        exit 1
    fi
fi

# Warning if compiler is not gfortran or gcc-based
if [[ "$FC" != *gfortran* && "$FC" != *gcc* ]]; then
    echo -e "${YELLOW}Warning: The code has been tested mainly with 'gfortran'.${NC}"
    echo -e "${YELLOW}         Using '$FC' may lead to unexpected issues.${NC}"
fi

# Handle build directory
if [[ -d "$build_dir" ]]; then
    read -p "The directory '$build_dir' already exists. Do you want to overwrite it? (y/n): " response
    response=$(echo "$response" | tr '[:upper:]' '[:lower:]')
    if [[ "$response" == "y" || "$response" == "yes" ]]; then
        rm -rf "$build_dir"
        echo -e "${GREEN}Directory '$build_dir' has been removed.${NC}"
    else
        echo -e "${RED}Operation cancelled.${NC}"
        exit 1
    fi
fi

mkdir "$build_dir"
cd "$build_dir" || exit 1


# Run CMake configuration
echo -e "${GREEN}Running CMake configuration...${NC}"
cmake .. -DCMAKE_BUILD_TYPE="$build_type" -DENABLE_OMP="$enable_omp"

if [[ $? -ne 0 ]]; then
    echo -e "${RED}CMake configuration failed.${NC}"
    exit 1
fi

echo -e "${GREEN}Configuration completed successfully!${NC}"
echo ""
echo "Summary:"
echo -e "  Build directory: ${YELLOW}${build_dir}${NC}"
echo -e "  Build type: ${YELLOW}${build_type}${NC}"
echo -e "  Fortran compiler: ${YELLOW}${FC}${NC}"
echo -e "  OpenMP enabled: ${YELLOW}${enable_omp}${NC}"
echo ""
echo "Next steps:"
echo "  cd $build_dir"
echo "  make"
echo ""
echo "To run tests:"
echo "  ctest"
echo ""

