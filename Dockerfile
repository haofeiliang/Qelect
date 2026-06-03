FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    cmake \
    libgmp3-dev \
    libntl-dev \
    git \
 && rm -rf /var/lib/apt/lists/*

# Proxy settings (pass via --build-arg, skipped if empty)
ARG HTTP_PROXY
ARG HTTPS_PROXY
RUN if [ -n "$HTTP_PROXY" ]; then git config --global http.proxy "$HTTP_PROXY"; fi \
 && if [ -n "$HTTPS_PROXY" ]; then git config --global https.proxy "$HTTPS_PROXY"; fi

# Set working directory
WORKDIR /home/ubuntu/qelect

# Install PALISADE library
RUN git clone --depth 1 -b v1.11.9 https://gitlab.com/palisade/palisade-release \
    && cd palisade-release && mkdir build && cd build \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/home/ubuntu/qelect/build \
    && make -j$(nproc) && make install

# Intel HEXL support (disable with --build-arg USE_INTEL_HEXL=OFF)
ARG USE_INTEL_HEXL=ON

# Install SEAL library
RUN git clone --depth 1 https://github.com/wyunhao/SEAL \
    && cd SEAL && cmake -S . -B build \
        -DCMAKE_INSTALL_PREFIX=/home/ubuntu/qelect/build \
        -DSEAL_USE_INTEL_HEXL=$USE_INTEL_HEXL \
    && cmake --build build -- -j$(nproc) && cmake --install build

# Copy the qelect repository
COPY . /home/ubuntu/qelect

# Build the qelect project
RUN mkdir -p /home/ubuntu/qelect/build /home/ubuntu/qelect/data/perm \
    && cd /home/ubuntu/qelect/build \
    && cmake .. -DCMAKE_PREFIX_PATH=/home/ubuntu/qelect/build \
    && make -j$(nproc)

# Set entrypoint
ENTRYPOINT ["/home/ubuntu/qelect/build/mps"]
