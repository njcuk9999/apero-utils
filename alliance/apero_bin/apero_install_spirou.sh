#!/bin/bash

# Set instrument
export APERO_INSTRUMENT="spirou"

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh
# pass through to global instrument script
source $APERO_BIN_PATH/apero_install_instrument.sh "$@"

# SPIROU specific post-install steps can be added below