# Start with a CentOS 7 base image
FROM centos:7

# Install necessary tools, libraries, and development tools
RUN yum install -y epel-release && \
    # Install necessary tools, libraries, and development tools
    yum install -y \
    atlas atlas-devel \
    gcc \
    gcc-c++ \
    python2-pip \
    python2-devel \
    tkinter python-tk \
    libhdf5-devel \
    netcdf-devel \
    lapack-devel \
    lapacke-devel \
    openblas openblas-devel \
    wget \
    libXext \
    mesa-libGL \
    libSM-devel && \
    yum groupinstall -y "Development Tools" && \
    # Clean up yum cache
    yum clean all && \
    rm -rf /var/cache/yum

# Set a working directory for the application
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Upgrade pip and setuptools
RUN pip install pip==20.3.4 && \
    pip install --upgrade 'setuptools>=18.0'

# Set environment variables for OpenBLAS
ENV LD_LIBRARY_PATH /usr/lib64/openblas:/usr/lib64/:$LD_LIBRARY_PATH
ENV LIBRARY_PATH /usr/lib64/openblas:$LIBRARY_PATH
ENV C_INCLUDE_PATH /usr/include/openblas:$C_INCLUDE_PATH

# Install the required Python packages
RUN pip install numpy scipy configparser cftime==1.0.4 pillow vtk matplotlib opencv-python==3.4.2.17 pytest Cython>=0.19 netCDF4==1.5.4

# Clean build artifacts to ensure a clean build and Install the package using setup.py
RUN python2.7 setup.py clean && \
    pip install .