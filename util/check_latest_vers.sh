#!/bin/sh

CUR_VERS=$(grep "VERSION" src/include/pihm.h |awk '{print $3}'|tr -d '"')
LATEST_VERS=$(curl --silent "https://api.github.com/repos/PSUmodeling/MM-PIHM/releases/latest" | grep '"tag_name":' | sed -E 's/.*"([^"]+)".*/\1/')
LATEST_VERS=${LATEST_VERS#"v"}

if [ "$(printf '%s\n' $CUR_VERS $LATEST_VERS | sort -V | head -n 1)" != "$LATEST_VERS" ]; then
    echo
    echo "MM-PIHM has been updated to v${LATEST_VERS}, and you are using v${CUR_VERS}."
    echo "You can download the latest version at https://github.com/PSUmodeling/MM-PIHM/releases/tag/v${LATEST_VERS}"
    echo
fi
