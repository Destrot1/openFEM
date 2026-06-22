#!/bin/bash

# Start the container if it is not already running
if [ ! "$(docker ps -q -f name=openFEM_dev)" ]; then
    echo "Starting container..."
    docker compose up -d
fi

echo "GUI viewer (gmsh 'view' command): http://localhost:6080/vnc.html"

# Enter the container, build, and run
docker exec -it openFEM_dev bash -c "
    cd /workspace/build
    cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON > /dev/null 2>&1
    make -j\$(nproc)
    if [ \$? -eq 0 ]; then
        echo 'Build OK'
        ./openFEM
    else
        echo 'Build failed'
    fi
"