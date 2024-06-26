# Use an official Python 2.7 runtime as a parent image
FROM python:2.7

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libgl1-mesa-glx \
    libhdf5-serial-dev \
    libnetcdf-dev \
    liblapack-dev \
    liblapacke-dev \
    libblas-dev \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables for HDF5 and NetCDF
ENV HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial
ENV NETCDF4_DIR=/usr/include

# Create and set a directory for the application
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Upgrade pip
RUN pip install --upgrade pip

# Install the required packages
RUN pip install -r requirements.txt

# Install the package using setup.py
RUN pip install .

RUN git clone https://github.com/xiyuanbao/openpiv-python-gpu.git
WORKDIR /app/openpiv-python-gpu
RUN python setup.py build_ext -i
RUN python setup.py install

WORKDIR /app
# Clone the pyflow repository
RUN git clone https://github.com/pathak22/pyflow.git

# Navigate to the cloned directory and build the extension
WORKDIR /app/pyflow
RUN python setup.py build_ext -i

WORKDIR /app