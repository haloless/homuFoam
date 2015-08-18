#!/bin/sh
# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

#
cleanCase
rm -rf logs > /dev/null 2>&1
rm -f *.eps > /dev/null 2>&1


