#!/bin/sh

CYCLES_DIR=$1
REQ_VERS=$2

# Check Cycles version
CUR_VERS=$(grep "VERSION" $CYCLES_DIR/include/cycles.h |awk '{print $3}'|tr -d '"')
PIHM_VERS=$(grep "VERSION" src/include/pihm.h |awk '{print $3}'|tr -d '"')

# Check if Cycles version is compatible
if [ "$(printf '%s\n' $CUR_VERS $REQ_VERS | sort -V | head -n 1)" != "$REQ_VERS" ]; then
    echo "MM-PIHM v$PIHM_VERS requires Cycles v$REQ_VERS or above."
    echo "Currently Cycles v$CUR_VERS is installed."
    exit 1
fi
