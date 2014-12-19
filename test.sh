#!/bin/bash

fulldir=`pwd`

function testit {
t=$1
fulldir=$2
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
  echo "!  $t ... FAIL";
  cat $t.tmp;
else

  sedthing="s#$fulldir#BURNMAN#g"
  sed -i'' -e $sedthing $t.tmp #remove the absolute path from warnings
  sed -i'' -e "s/.py:[0-9]*:/.py:/g" $t.tmp #remove line numbers
  sed -i'' -e '/UserWarning: findfont: Could not match/d' $t.tmp #remove font warning crap
  sed -i'' -e '/UserWarning: findfont: Font family/d' $t.tmp #remove font warning crap
  sed -i'' -e '/tight_layout : falling back to Agg renderer/d' $t.tmp #remove font warning crap
  sed -i'' -e '/cannot be converted with the encoding. Glyph may be wrong/d' $t.tmp #remove font warning crap
  sed -i'' -e '/time old .* time new/d' $t.tmp #remove timing from tests/debye.py

  (numdiff -r 1e-5 -s ' \t\n[],' -a 1e-5 -q $t.tmp $fulldir/misc/ref/$t.out >/dev/null && rm $t.tmp && echo "  $t ... ok"
  ) || {
  echo "!  $t ... FAIL";
  echo "Check: `readlink -f $t.tmp` `readlink -f $fulldir/misc/ref/$t.out`";
  numdiff -r 1e-5 -s ' \t\n[],' -a 1e-5 $t.tmp $fulldir/misc/ref/$t.out | head
  }

fi
}




echo "*** running test suite..."

# check for tabs in code:
for f in `find . -name \*.py`
do
    grep -P "\t" -q $f && echo "ERROR: tabs found in '$f'" && exit 1
done

cd tests
python tests.py || exit 1
cd ..


cd misc
echo "gen_doc..."
python gen_doc.py >/dev/null || exit 1

cd benchmarks
for test in `ls *.py`
do
    testit $test $fulldir
done
cd ..
cd ..




echo "checking examples/ ..."
cd examples
for test in `ls example*.py`
do
    [ $test == "example_inv_murakami.py" ] && echo "  *** skipping $test !" && continue
    [ $test == "example_premite_isothermal.py" ] && echo "  *** skipping $test !" && continue
    [ $test == "example_compare_enstpyro.py" ] && echo "  *** skipping $test !" && continue
    [ $test == "example_partition_coef.py" ] && echo "  *** skipping $test !" && continue

    testit $test $fulldir
done
cd ..


echo "checking misc/ ..."
cd misc
for test in `ls paper*.py`
do
    [ $test == "paper_opt_pv_old.py" ] && echo "  *** skipping $test !" && continue

    testit $test $fulldir
done

testit table.py $fulldir
cd ..

echo "   done"



echo ""
echo "*** tests done"
