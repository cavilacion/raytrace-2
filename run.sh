#!/bin/bash
echo "Configuring project..."
cd build
cmake ..
echo "Building project..."
make
echo "Running project..."
./ray ../Scenes/scene01-lights-shadows.json
echo "Finished."
