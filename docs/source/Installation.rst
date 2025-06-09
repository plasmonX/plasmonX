.. _Installation:

Installation
============

Choose your platform:

- `Install on Linux <#install-on-linux>`__
- `Install on macOS <#install-on-macos>`__


Install on Linux
----------------

1. Clone the repository:

   .. code-block:: bash

      git clone --recursive https://github.com/your_project.git

2. Install python:

   2.1. If you are on a local system, the suggested option is to install python and all the requirements by:

        .. code-block:: bash

           sudo apt update
           sudo apt install python3 python3-pip

        .. code-block:: bash

           pip install -r requirements.txt

   2.2. Alternatively, use conda:

        .. code-block:: bash

           conda env create -f python/requirements.yaml
           conda activate plasmonX-env

        If you don't have conda, you can install Miniconda with:

        .. code-block:: bash

           wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
           bash Miniconda3-latest-Linux-x86_64.sh
           source ~/.bashrc

3. (Suggested) Install MKL for Lapack/BLAS

4. Install CMake [minumum version required: 3.15]

   .. code-block:: bash
   
      sudo apt get install cmake

5. Run the setup.sh [type -h for options]

   .. code-block:: bash
   
      ./setup.sh -omp -b BUILD_DIR

6. Compile the code

   .. code-block:: bash
   
      cd BUILD_DIR
      make -j

7. Run the tests [Please, avoid the -j option]

   .. code-block:: bash
   
      ctest

Install on macOS
----------------

1. Clone the repository:

   .. code-block:: bash

      git clone --recursive https://github.com/your_project.git

2. Install python:

   2.1. If you are on a local system, the suggested option is to install python and all the requirements by:

        .. code-block:: bash

           brew install python

        .. code-block:: bash

           pip3 install -r requirements.txt

   2.2. Alternatively, use conda:

        .. code-block:: bash

           conda env create -f python/requirements.yaml
           conda activate plasmonX-env

        If you don't have conda, you can install Miniconda with:

        .. code-block:: bash

           curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
           bash Miniconda3-latest-MacOSX-x86_64.sh
           source ~/.bash_profile

3. (Suggested) Install MKL for Lapack/BLAS

4. Install CMake [minimum version required: 3.15]

   .. code-block:: bash

      brew install cmake

5. Run the setup.sh [type -h for options]

   .. code-block:: bash

      ./setup.sh -omp -b BUILD_DIR

6. Compile the code

   .. code-block:: bash

      cd BUILD_DIR
      make -j

7. Run the tests [Please, avoid the -j option]

   .. code-block:: bash

      ctest
