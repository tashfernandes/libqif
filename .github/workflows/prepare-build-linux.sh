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

if [ -n "$PYTHON_LIBRARY" ]; then
    EXTRA_PY_ARGS="-DPython3_LIBRARY=$PYTHON_LIBRARY"
else
    echo "No libpython found, continuing without it"
    EXTRA_PY_ARGS=""
fi

mkdir build
cd build
echo "About to run... "
echo "cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DPython3_EXECUTABLE=$PYTHON \
      -DPython_FIND_STRATEGY=LOCATION \
      -DPython3_INCLUDE_DIR=$PYTHON_INCLUDE_DIR \
      -DPython3_FIND_DEBUG=ON \
      $EXTRA_PY_ARGS .."
#make install -j 2
