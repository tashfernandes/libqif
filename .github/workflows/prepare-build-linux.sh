#!/bin/bash

set -x		# echo commands
set -e		# stop on any error

PYTHON=$(which python3.12)
echo $PYTHON
cmake --version

PYTHON_INCLUDE_DIR=$($PYTHON -c "from sysconfig import get_paths; print(get_paths()['include'])")
PYTHON_LIBRARY=$(find $(dirname $PYTHON)/../lib -name "libpython3.12*.so" | head -n 1)
echo $PYTHON_INCLUDE_DIR
echo $PYTHON_LIBRARY

export PATH=/opt/python/cp312-cp312/bin:$PATH
export CMAKE_POLICY_VERSION_MINIMUM=3.5

cat /proc/cpuinfo

if [ -n "$PYTHON_LIBRARY" ]; then
    EXTRA_PY_ARGS="-DPYTHON_LIBRARY=$PYTHON_LIBRARY"
else
    echo "No libpython found, continuing without it"
    EXTRA_PY_ARGS=""
fi

mkdir build
cd build

cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DPYBIND11_FINDPYTHON=ON \
      -DPYTHON_EXECUTABLE=$PYTHON \
      -DPYTHON_ROOT_DIR=/opt/python/cp312-cp312 \
      -DPYTHON_INCLUDE_DIR=$PYTHON_INCLUDE_DIR \
      $EXTRA_PY_ARGS ..
make install -j 2
