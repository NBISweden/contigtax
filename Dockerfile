FROM continuumio/miniconda3:4.8.2

LABEL maintainer="John Sundh" email=john.sundh@nbis.se
LABEL description="Docker image for python package contigtax"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /work

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && apt-get clean

# Add environment file
COPY environment.yaml .

# Install environment into base
RUN conda env update -n base -f environment.yaml && conda clean -a

# Install tango
COPY contigtax contigtax
COPY requirements.txt .
COPY setup.py .
COPY README.md .
RUN python -m pip install .

# Set entrypoint
ENTRYPOINT ["contigtax"]