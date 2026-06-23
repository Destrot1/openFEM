FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    build-essential \
    cmake \
    git \
    curl \
    zip \
    unzip \
    tar \
    pkg-config \
    ninja-build \
    libgmsh-dev \
    gmsh \
    libocct-foundation-dev \
    libocct-modeling-algorithms-dev \
    libocct-modeling-data-dev \
    libocct-data-exchange-dev \
    libocct-visualization-dev \
    libeigen3-dev \
    nlohmann-json3-dev \
    libfmt-dev \
    ttyd \
    libgl1 \
    libxext6 \
    libxrender1 \
    libxft2 \
    libfontconfig1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

COPY docker/start-gui.sh /usr/local/bin/start-gui.sh
RUN chmod +x /usr/local/bin/start-gui.sh


ENTRYPOINT ["/usr/local/bin/start-gui.sh"]
CMD ["/bin/bash"]