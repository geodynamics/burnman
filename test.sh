#!/bin/bash

function testit {
t=$1
echo "*** testing $t..."
(python <<EOF
import matplotlib as m
m.use('Cairo')
VERBOSE=1
execfile('$t')
EOF
)  >$t.tmp 2>&1 || { echo "test $t failed"; exit 1; } 
echo "diff $t.tmp misc/ref/$t.out"
diff $t.tmp misc/ref/$t.out || { echo "test $t failed"; exit 1; } 
rm $t.tmp
}




echo "*** running test suite..."
cd tests
python tests.py || exit 1
cd ..

python burnman/composition.py || exit 1

cd misc
python gen_doc.py >/dev/null || exit 1
cd ..

testit "burnman/benchmark.py"
echo "   done"


echo ""

t="example_composition.py"
testit $t
diff example_composition.png misc/ref/example_composition.png || { echo "test $t failed"; exit 1; } 
echo "   done"

testit "example_compare_enstpyro.py"
echo "   done"

testit "example_seismic.py"
echo "   done"



testit "example_seismic.py"
echo "   done"
testit "example_compare_two_models.py"
echo "   done"
testit "example_spintransition.py"
echo "   done"
testit "example_composition.py"
echo "   done"
testit "example_user_input_material.py"
echo "   done"
testit "example_geotherms.py"
echo "   done"

testit "example_woutput.py"
echo "   done"

testit "example_optimize_pv.py"
echo "   done"
testit "table.py"
echo "   done"



echo ""
echo "*** tests done"
