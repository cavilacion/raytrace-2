#!/bin/bash
echo "Configuring project..."
cd build
cmake ..
echo "Building project..."
make
echo "Running project..."
./ray ../Scenes/example.json
echo "Finished."
