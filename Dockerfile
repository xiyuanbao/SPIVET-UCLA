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
    libhdf5-devel \
    netcdf-devel \
    lapack-devel \
    lapacke-devel \
    wget && \
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

# Install the required Python packages
RUN pip install numpy scipy configparser cftime==1.0.4 pillow vtk matplotlib opencv-python==3.4.2.17 pytest Cython>=0.19 netCDF4==1.5.4

# Install the package using setup.py
RUN pip install .

# Define the command to run on container start
CMD ["pytest", "tests/."]