#!/bin/bash

fulldir=`pwd`

# use numdiff but fall back to diff if not found:
which numdiff >/dev/null
if [ $? -eq 0 ]
then
  diffcmd="numdiff -r 2e-5 -s ' \t\n[],'"
else
  diffcmd="diff"
  echo "WARNING: numdiff not found, please install! Falling back to diff."
fi

if [ -z $PYTHON ]
then
  PYTHON=python
fi

$PYTHON --version

function testit {
t=$1

if [ "$t" == "__init__.py" ]
then
return
fi

fulldir=$2
#echo "*** testing $t ..."
($PYTHON <<EOF
import matplotlib as m
m.use('Template')
VERBOSE=1
RUNNING_TESTS=1
with open('$t') as f:
    CODE = compile(f.read(), '$t', 'exec')
    exec(CODE)
EOF
) >$t.tmp 2>$t.tmp.error
ret=$?
cat $t.tmp.error >>$t.tmp
rm -f $t.tmp.error
if [ "$ret" -ne 0 ]
then
  echo "!  $t ... FAIL";
  cat $t.tmp;
else

  sedthing="s#$fulldir#BURNMAN#g"
  sed -i.bak -e $sedthing $t.tmp #remove the absolute path from warnings
  sed -i.bak -e "s/.py:[0-9]*:/.py:/g" $t.tmp #remove line numbers
  sed -i.bak -e 's/<stdin>:[0-9]*: //g' $t.tmp #remove line numbers
  sed -i.bak -e '/UserWarning: findfont: Could not match/d' $t.tmp #remove font warning crap
  sed -i.bak -e '/UserWarning: findfont: Font family/d' $t.tmp #remove font warning crap
  sed -i.bak -e '/tight_layout : falling back to Agg renderer/d' $t.tmp #remove font warning crap
  sed -i.bak -e '/cannot be converted with the encoding. Glyph may be wrong/d' $t.tmp #remove font warning crap
  sed -i.bak -e '/time old .* time new/d' $t.tmp #remove timing from tests/debye.py
  sed -i.bak -e '/  relative central core pressure error.*/d' $t.tmp #remove residuals from examples/example_build_planet

  rm -f $t.tmp.bak

  (eval $diffcmd $t.tmp $fulldir/misc/ref/$t.out >/dev/null && rm $t.tmp && echo "  $t ... ok"
  ) || {
  echo "!  $t ... FAIL";
  echo "Check: `pwd`/$t.tmp $fulldir/misc/ref/$t.out";
  eval $diffcmd $t.tmp $fulldir/misc/ref/$t.out | head
  }

fi
}




echo "*** running test suite..."

# check for tabs in code:
for f in `find . -name \*.py | grep -v ipython/`
do
    
    grep $'\t' -q $f && \
	echo "ERROR: tabs found in '$f':" && \
	grep -n $'\t' $f && exit 0
done


cd tests
$PYTHON tests.py || (echo "ERROR: unittests failed"; exit 1) || exit 0
cd ..


cd misc
echo "gen_doc..."
$PYTHON gen_doc.py >/dev/null || exit 0

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
    [ $test == "example_seismic_travel_times.py" ] && echo "  *** skipping $test !" && continue
    testit $test $fulldir
done
cd ..


echo "checking misc/ ..."
cd misc
for test in `ls *.py`
do
    [ $test == "gen_doc.py" ] && echo "  *** skipping $test !" && continue
    [ $test == "table.py" ] && echo "  *** skipping $test !" && continue
    [ $test == "create_burnman_readable_perplex_table.py" ] && echo "  *** skipping $test !" && continue
    [ $test == "helper_solid_solution.py" ] && echo "  *** skipping $test !" && continue

    testit $test $fulldir
done
testit table.py $fulldir
cd ..

echo "checking contrib/CHRU2014 ..."
cd contrib/CHRU2014
for test in `ls *.py`
do
    testit $test $fulldir
done
cd ../..


echo "checking contrib/tutorial/ ..."
cd contrib/tutorial/
for test in `ls step*.py`
do
testit $test $fulldir
done
cd ../..


echo "checking contrib/CHRU2014/ ..."
cd contrib/CHRU2014
for test in `ls *.py`
do
testit $test $fulldir
done
cd ../..


echo "   done"

echo ""
echo "*** tests done"
