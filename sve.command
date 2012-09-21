#!/bin/bash
# Directory of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Default exeuction
SVE='/usr/bin/env python '$DIR'/interactive.py'

if [ ! -f interactive.py ]; then
    platform=`uname`
    if [[ "$platform" == 'Linux' && -f $DIR/sve.linux2/sve.linux2 ]]; then
        SVE=$DIR'/sve.linux2/sve.linux2'
    elif [[ "$platform" == 'Darwin' && -f $DIR/sve.darwin/sve.darwin ]]; then
        SVE=$DIR'/sve.darwin/sve.darwin'
    fi
fi

$SVE "$@"
