#!/bin/bash


name="burnman-v0.3"

rm -rf $name

svn export https://burnman.googlecode.com/svn/trunk $name
rm -rf $name/misc/ $name/todo.txt $name/release.sh $name/test.sh $name/reproduce_*.py $name/play_withslb.py



zip -v $name.zip $name
