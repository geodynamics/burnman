#!/bin/bash

fulldir=`pwd`

function testit {
t=$1
#echo "*** testing $t ..."
(python <<EOF
import matplotlib as m
m.use('Template')
VERBOSE=1
RUNNING_TESTS=1
execfile('$t')
EOF
) >$t.tmp 2>&1
ret=$?
if [ "$ret" -ne 0 ]
then
  echo "$t ... FAIL";
  cat $t.tmp;
else

  sedthing="s#$fulldir#BURNMAN#g"
  sed -i'' $sedthing $t.tmp #remove the absolute path from warnings
  sed -i'' "s/.py:[0-9]*:/.py:/g" $t.tmp #remove line numbers
  sed -i'' '/UserWarning: findfont: Could not match/d' $t.tmp #remove font warning crap
  sed -i'' '/UserWarning: findfont: Font family/d' $t.tmp #remove font warning crap
  sed -i'' '/tight_layout : falling back to Agg renderer/d' $t.tmp #remove font warning crap
  sed -i'' '/cannot be converted with the encoding. Glyph may be wrong/d' $t.tmp #remove font warning crap
  sed -i'' '/time old .* time new/d' $t.tmp #remove timing from tests/debye.py

  (numdiff -r 1e-6 -s ' \t\n[]' -a 1e-6 -q $t.tmp ../misc/ref/$t.out >/dev/null && rm $t.tmp && echo "  $t ... ok"
  ) || {
  echo "$t ... FAIL";
  echo "Check: `readlink -f $t.tmp` `readlink -f ../misc/ref/$t.out`";
  numdiff -r 1e-6 -s ' \t\n[]' -a 1e-6 $t.tmp ../misc/ref/$t.out | head
  }

fi
}




echo "*** running test suite..."
cd tests
python tests.py || exit 1
testit "benchmark.py"
testit "debye.py"
cd ..

cd misc
echo "gen_doc..."
python gen_doc.py >/dev/null || exit 1
cd ..




echo "checking examples/ ..."
cd examples
for test in `ls example*.py`
do
    [ $test == "example_inv_murakami.py" ] && echo "*** skipping $test !" && continue
    [ $test == "example_premite_isothermal.py" ] && echo "*** skipping $test !" && continue
    [ $test == "example_compare_enstpyro.py" ] && echo "*** skipping $test !" && continue
    [ $test == "example_partition_coef.py" ] && echo "*** skipping $test !" && continue

    testit $test
done
cd ..


echo "checking misc/ ..."
cd misc
for test in `ls paper*.py`
do
    [ $test == "paper_opt_pv_old.py" ] && echo "*** skipping $test !" && continue

    testit $test
done

testit table.py
cd ..

echo "   done"



echo ""
echo "*** tests done"
