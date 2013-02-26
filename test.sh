#!/bin/bash

tests="example_compare_enstpyro.py example_seismic.py example_compare_two_models.py example_spintransition.py example_composition.py example_user_input_material.py example_geotherms.py example_woutput.py example_optimize_pv.py"
for t in $tests; do
    echo "*** testing '$t' ***"
    python $t || { echo "test $t failed"; exit 1; } 
done

echo ""
echo "*** tests done"
