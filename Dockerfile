# syntax=docker/dockerfile:1

# Use a base image with Miniconda
FROM continuumio/miniconda3:latest

# Set the working directory as base
WORKDIR /app

# Install app dependencies and git necessary for submodules when using info from ecoli_fbi github repository
RUN apt-get update && apt-get install -y git  && \
    /opt/conda/bin/conda install -c conda-forge pytest 

RUN pip install bifrostlib

# Copy the entire repository into the container
COPY . .

# Copy the install.sh and environment.yml into the container
COPY install.sh .
COPY environment.yml .

# Ensure install.sh is executable
RUN chmod +x install.sh

# Initialize conda for bash shell
RUN /opt/conda/bin/conda init bash

# Install the tool using the install script
RUN bash install.sh -i LOCAL

# Set environment variable for Conda environment name
ARG CONDA_ENV_NAME=""
ENV CONDA_ENV_NAME="${CONDA_ENV_NAME}"

# Set the default command to activate the Conda environment and run pytest
CMD ["/bin/bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate ${CONDA_ENV_NAME} && python -m pytest"]