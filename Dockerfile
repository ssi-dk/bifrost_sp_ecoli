# syntax=docker/dockerfile:1

# Use a base image with Miniconda
FROM continuumio/miniconda3:latest

# Set the working directory as base
WORKDIR /app

# Install app dependencies and git necessary for submodules when using info from ecoli_fbi github repository
RUN apt-get update && apt-get install -y git

# Copy the entire repository into the container
COPY . .

# Initialize and update submodules
RUN git submodule init && git submodule update

# Copy the install.sh and environment.yml into the container
COPY environment.yml .
COPY install.sh .

# Ensure install.sh is executable
RUN chmod +x install.sh

# Install dependencies and create the conda environment
RUN conda env create -f environment.yml

# Set the default shell to conda environment
SHELL ["conda", "run", "-n", "bifrost_sp_ecoli_env", "/bin/bash", "-c"]

# Install the tool using the install script
RUN bash install.sh -i LOCAL

# Set environment variables
ENV BIFROST_INSTALL_DIR='/app'
# Use ARG for database key and set at runtime
ARG BIFROST_DB_KEY=''

# Set the default command to run the Python module
CMD ["bash"]

# CMD ["bash", "-c", "conda activate bifrost_sp_ecoli_env && pytest"] consider changing the workflow to accomodate one way or the other, by solely using pytest
