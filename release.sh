#!/bin/bash


name="burnman-v0.3"

rm -rf $name

svn export https://burnman.googlecode.com/svn/trunk $name
rm -rf $name/misc/ $name/todo.txt $name/release.sh $name/test.sh play_withslb.py $name/reproduce_*.py $name/play_withslb.py



zip -v $name.zip *py readme.txt burnman/*py data/*.txt
