FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
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
    xvfb \
    x11vnc \
    novnc \
    websockify \
    ttyd \
    libgl1 \
    libxext6 \
    libxrender1 \
    libxft2 \
    libfontconfig1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Starts a virtual display + VNC + noVNC web bridge so GUI windows
# opened inside the container (e.g. gmsh's viewer) can be seen from a
# browser on the host at http://localhost:6080 -- works the same on
# any host OS, no X server needed on the host.
COPY docker/start-gui.sh /usr/local/bin/start-gui.sh
RUN chmod +x /usr/local/bin/start-gui.sh

ENV DISPLAY=:99

ENTRYPOINT ["/usr/local/bin/start-gui.sh"]
CMD ["/bin/bash"]