# Start with a CentOS 7 base image
FROM centos:7

# Install EPEL Repository for additional packages
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
    rm -rf /var/cache/yum \

# Create and set a directory for the application
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Upgrade pip, setuptools and install required packages
RUN pip install pip==20.3.4 && \
    pip install --upgrade 'setuptools>=18.0' && \
    pip install \
    numpy \
    scipy \
    configparser \
    cftime==1.0.4 \
    pillow \
    vtk \
    matplotlib \
    opencv-python==3.4.2.17 \
    pytest \
    Cython>=0.19 && \
    pip download netCDF4==1.5.4 --no-deps && \
    tar -xzvf netCDF4-1.5.4.tar.gz && \
    cd netCDF4-1.5.4 && \
    python setup.py egg_info

# Install the package using setup.py
RUN pip install .

# Define the command to run on container start
CMD ["pytest", "tests/."]