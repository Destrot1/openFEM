#!/bin/bash

# Create a shortcut so STEP files can be loaded as /Desktop/foo.step instead of
# /host_home/Desktop/foo.step
ln -sfn /host_home/Desktop /Desktop

# Delete meshes/results of the previous session
rm -f /workspace/web/scene.json

# Serve the viewer page on
# http://localhost:8081/ from any browser.
( cd /workspace/web && python3 -m http.server 8081 ) &

# Run cmake to compile autonomously the compilation of c++ scripts
mkdir -p /workspace/build
( cd /workspace/build && \
    { [ -f CMakeCache.txt ] || cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON > /dev/null 2>&1 ; } && \
    cmake --build . -j$(nproc) 
)
# A real interactive shell in the browser (ttyd + xterm.js), embedded
# directly in the web viewer page (but served on port 7681). Visual 
# integration done by <iframe> from the 8081 page (<iframe src="http://localhost:7681/"> )
( cd /workspace && ttyd --writable \
    -t fontFamily='"JetBrains Mono, monospace"' \
    -t fontSize=13 \
    -t 'theme={"background":"#141a22","foreground":"#eef1f6","cursor":"#5b8def","selectionBackground":"rgba(91,141,239,0.3)"}' \
    -p 7681 bash -c 'cd /workspace/build && ./openFEM; exec bash' ) &

exec "$@"
