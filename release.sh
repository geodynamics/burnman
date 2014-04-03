#!/bin/bash


name="burnman-v0.6"

rm -rf $name $name.zip

svn export https://burnman.googlecode.com/svn/trunk $name
rm -rf $name/homepage  $name/todo.txt $name/release.sh $name/test.sh $name/reproduce_*.py $name/play_withslb.py $name/example_inv_big_pv.py $name/example_inv_murakami.py $name/example_optimize_slb2011.py $name/example_premite_isothermal.py $name/table_latex.py $name/sigfig_table.py $name/misc/ref $name/tests $name/murakami_book_chapter.py $name/input_figures/matas* $name/data $name/input_figures/murakami* $name/burnman/benchmark.py
#$name/misc/

zip -rv $name.zip $name
