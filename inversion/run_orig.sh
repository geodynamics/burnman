#!/bin/bash

for i in "1" "2" "3" "4";
do
  f=orig_$i
  if [ -e $f ]
  then
    python setup_isochemical_newmisfit_moduli_realdata.py continue $f &
  else
    python setup_isochemical_newmisfit_moduli_realdata.py run $f &
  fi
done

wait
