#!/bin/bash


name="burnman-v0.3"

rm -rf $name

svn export https://burnman.googlecode.com/svn/trunk $name
rm -rf $name/misc/ $name/todo.txt $name/release.sh $name/test.sh $name/reproduce_*.py $name/play_withslb.py $name/example_inv_big_pv.py $name/example_inv_murakami.py $name/example_optimize_slb2011.py $name/example_premite_isothermal.py

zip -rv $name.zip $name
