#!/bin/bash

echo -e "This script will install dependencies and may modify or break system packages."
echo "If you prefer, you can manually install dependencies by running:"
echo "    python3 -m venv geom_venv && source geom_venv/bin/activate"
echo "    pip install --upgrade pip && pip install -r requirements.txt"
echo

echo -n "Do you want to proceed with the installation? (Y/n): "
read -r response
if [[ "$response" =~ ^([nN][oO]?|[nN])$ ]]; then
    echo "Installation aborted."
    exit 0
fi

echo "Setting up GEOM project..."

# Ensure the script stops on any error
set -e

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Detect OS type
OS_TYPE="$(uname -s)"
if [[ "$OS_TYPE" == "Darwin" ]]; then
    OS="macOS"
elif [[ "$OS_TYPE" == "Linux" ]]; then
    if grep -qi microsoft /proc/version; then
        OS="WSL"
    elif command_exists apt; then
        OS="Debian"
    elif command_exists dnf; then
        OS="Fedora"
    else
        OS="Linux-Other"
    fi
elif [[ "$OS_TYPE" =~ "MINGW" || "$OS_TYPE" =~ "CYGWIN" || "$OS_TYPE" =~ "MSYS" ]]; then
    OS="Windows"
else
    echo "Unsupported OS: $OS_TYPE"
    exit 1
fi

echo "Detected OS: $OS"

# Install Python and Pip
if [[ "$OS" == "macOS" ]]; then
    echo "Using Homebrew to install Python..."
    if ! command_exists brew; then
        echo "Error: Homebrew is not installed. Install it from https://brew.sh/"
        exit 1
    fi
    if ! command_exists python3; then
        brew install python
    fi
    python3 -m pip install --upgrade pip setuptools wheel --break-system-packages

elif [[ "$OS" == "Debian" ]]; then
    echo "Using APT to install Python..."
    sudo apt update
    sudo apt install -y python3 python3-pip

elif [[ "$OS" == "Fedora" ]]; then
    echo "Using DNF to install Python..."
    sudo dnf install -y python3 python3-pip

elif [[ "$OS" == "WSL" ]]; then
    echo "Detected WSL. Using APT (or alternative) to install Python..."
    if command_exists apt; then
        sudo apt update
        sudo apt install -y python3 python3-pip
    elif command_exists dnf; then
        sudo dnf install -y python3 python3-pip
    else
        echo "No known package manager found. Please install Python manually."
        exit 1
    fi

elif [[ "$OS" == "Windows" ]]; then
    echo "Using Windows package manager to install Python..."
    if ! command_exists python; then
        echo "Installing Python using winget..."
        winget install -e --id Python.Python
    fi
    python -m pip install --upgrade pip

else
    echo "Unsupported Linux distribution. Please install Python manually."
    exit 1
fi

# Install project dependencies
echo "Installing Python dependencies..."
pip3 install --upgrade pip setuptools wheel --break-system-packages
pip3 install -r requirements.txt --break-system-packages

# Determine the correct shell configuration file
if [[ "$OS" == "Windows" ]]; then
    SHELL_RC="$HOME/Documents/PowerShell/Microsoft.PowerShell_profile.ps1"
elif [[ "$SHELL" == "/bin/zsh" ]]; then
    SHELL_RC="$HOME/.zshrc"
elif [[ "$SHELL" == "/bin/bash" ]]; then
    SHELL_RC="$HOME/.bashrc"
elif [[ "$SHELL" == "/bin/fish" ]]; then
    SHELL_RC="$HOME/.config/fish/config.fish"
else
    SHELL_RC="$HOME/.profile"
fi

# Define the geom_load function
GEOM_LOAD_FUNCTION="function geom_load {
    export PYTHONPATH=\"$PWD:\$PYTHONPATH\"
    alias geom=\"python3 -m geom\"
}"

# Add geom_load function if not already present
if ! grep -q "function geom_load" "$SHELL_RC"; then
    echo "Adding geom_load function to $SHELL_RC..."
    echo -e "\n$GEOM_LOAD_FUNCTION" >> "$SHELL_RC"
fi

# Run tests
if [ -f "./geom/tests/run_all_tests.sh" ]; then
    echo "Running tests..."
    chmod +x ./geom/tests/run_all_tests.sh
    ./geom/tests/run_all_tests.sh
else
    echo "No test script found at ./geom/tests/run_all_tests.sh. Skipping tests."
fi

# Final message
echo -e "\n\033[1;32m✔ Installation complete!\033[0m"
echo -e "\n\033[1;34m➡ Before using \`geom\`, you must first run:\033[0m"
echo -e "\033[1;33msource $SHELL_RC\033[0m   \033[0;37m# This sets up your environment\033[0m"
echo -e "\n\033[1;34m➡ After that, you can use:\033[0m"
echo -e "\033[1;33mgeom_load\033[0m   \033[0;37m# Loads the environment to use \`geom\` as a command\033[0m"
echo -e "\033[1;33mgeom -h\033[0m   \033[0;37m  # Shows available options\033[0m"
echo -e "\n\033[1;34m➡ Alternatively, you can always run:\033[0m"
echo -e "\033[1;33mpython3 -m geom -h\033[0m   \033[0;37m# Runs GEOM without loading the environment\033[0m\n"
