#!/bin/bash

echo "Installing dependencies..."

# Function to check and install pBLAT for the correct OS
install_pblat() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "Detected Linux. Installing Linux version of pBLAT..."
        wget https://anaconda.org/bioconda/pblat/2.5.1/download/linux-64/pblat-2.5.1-h84c94e8_2.tar.bz2
        conda install ./pblat-2.5.1-h84c94e8_2.tar.bz2
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo "Detected macOS. Installing macOS version of pBLAT..."
        wget https://anaconda.org/bioconda/pblat/2.5.1/download/osx-64/pblat-2.5.1-ha5265cf_2.tar.bz2
        conda install ./pblat-2.5.1-ha5265cf_2.tar.bz2
    else
        echo "Unsupported OS type: $OSTYPE"
        exit 1
    fi
}

# Check and install GCC
if ! command -v gcc &> /dev/null; then
    echo "GCC not found. Installing..."
    sudo apt-get install -y gcc
else
    echo "GCC is already installed."
fi

# Check and install R
if ! command -v R &> /dev/null; then
    echo "R not found. Installing..."
    sudo apt-get install -y r-base
else
    echo "R is already installed."
fi

# Check and install Python3 and pip
if ! command -v python3 &> /dev/null; then
    echo "Python3 not found. Installing..."
    sudo apt-get install -y python3 python3-pip
else
    echo "Python3 is already installed."
fi

# Install Python dependencies
if [ -f requirements.txt ]; then
    pip3 install -r requirements.txt
fi

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda not found. Installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    conda init
    source ~/.zshrc  # Use ~/.zshrc for macOS
else
    echo "Conda is already installed."
fi

# Install pBLAT for the correct OS
install_pblat

echo "Installation complete."
