#!/usr/bin/env bash

cdir=$(pwd) # full path to setup dir
param_dir=$cdir/paramFiles

for run in $(ls $cdir/../runs); do
  for f in $(ls $param_dir); do
    ln -sf $param_dir/$f $cdir/../runs/$run
  done
done
