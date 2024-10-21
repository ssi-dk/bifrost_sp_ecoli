# syntax=docker/dockerfile:1

# Use a base image with Miniconda
FROM continuumio/miniconda3:latest

# Set the working directory as base
WORKDIR /app

# Install app dependencies and git necessary for submodules when using info from ecoli_fbi github repository
RUN apt-get update && apt-get install -y git

# Copy the entire repository into the container
COPY . .

# Copy the install.sh and environment.yml into the container
COPY install.sh .
#COPY environment.yml .

# Ensure install.sh is executable
RUN chmod +x install.sh

# Initialize conda for bash shell
RUN /opt/conda/bin/conda init bash

# Install the tool using the install script
RUN bash install.sh -i LOCAL

# Set environment variables
ENV BIFROST_INSTALL_DIR='/app'
# Use ARG for database key and set at runtime
ARG BIFROST_DB_KEY=''

# Set the default command to run the Python module
CMD ["bash", "-c", "BIFROST_DB_KEY=\"$BIFROST_DB_KEY\" conda run -n $CONDA_ENV_NAME python cmd_test.py"]

