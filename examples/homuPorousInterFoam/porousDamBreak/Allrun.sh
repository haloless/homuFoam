#!/bin/sh
# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

application=`getApplication`

# generate fake for visualization
THISCASE=`basename $PWD`
touch $THISCASE.OpenFOAM

# copy the original state
cp 0/alpha.water.org 0/alpha.water
cp 0/varepsilon.org 0/varepsilon

# generate mesh
runApplication blockMesh
runApplication checkMesh

# set initial state
setFields

# run
runApplication $application


#runApplication snappyHexMesh -overwrite

#rm -f constant/{cellLevel,pointLevel}

#runApplication checkMesh -constant
#runApplication renumberMesh -overwrite

# Run
#runApplication $application

# Sample
#runApplication foamLog log.simpleFoam

# Postprocessing
#gnuplot res.gp



