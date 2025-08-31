#!/bin/bash

set -x		# echo commands
set -e		# stop on any error

PYTHON=$(which python3)
echo $PYTHON
cmake --version

PYTHON_INCLUDE_DIR=$($PYTHON -c "from sysconfig import get_paths; print(get_paths()['include'])")
PYTHON_LIBRARY=$(find $(dirname $PYTHON)/../lib -name "libpython3.12*.so" | head -n 1)
echo $PYTHON_INCLUDE_DIR
echo $PYTHON_LIBRARY

cat /proc/cpuinfo

mkdir build
cd build

cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DPython3_EXECUTABLE=$PYTHON \
      -DPython_FIND_STRATEGY=LOCATION \
      -DPython3_INCLUDE_DIR=$PYTHON_INCLUDE_DIR \
      -DPython3_FIND_DEBUG=ON \
      -DPython3_LIBRARY=$PYTHON_LIBRARY ..
#make install -j 2
