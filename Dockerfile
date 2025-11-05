# GPF Protein Design - Reproducible Environment
# Base: Ubuntu 22.04 (LTS)
FROM ubuntu:22.04

# Prevent interactive prompts during package install
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    git \
    build-essential \
    python3 \
    python3-pip \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Install FoldX (academic license required for commercial use)
RUN mkdir -p /opt/foldx && cd /opt/foldx && \
    wget http://foldxsuite.crg.eu/system/edd/5/FoldX5_Linux.tar.gz && \
    tar -xzf FoldX5_Linux.tar.gz && \
    rm FoldX5_Linux.tar.gz && \
    chmod +x foldx

# Add FoldX to PATH
ENV PATH="/opt/foldx:${PATH}"

# Set working directory
WORKDIR /app

# Copy requirements first (for layer caching)
COPY requirements.txt .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy entire repo
COPY . .

# Run unit tests during build (fail if tests fail)
RUN python3 -m pytest tests/ -v

# Default command
CMD ["python3", "-c", "print('GPF ready. Try: from gpf import gpf_transform_v4')"]
