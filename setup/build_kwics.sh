#!/usr/bin/env bash

base_dir="$(pwd)"
build_dir="$base_dir/../build"        # contains build folders for each theta
exec_dir="$base_dir/../exec"          # to put executables when done
mod_dir="$base_dir/../code"           # directory containing modified code
mitgcm_root="$HOME/MITgcm/MITgcm"     # MITgcm folder (from github)
buildopts="$HOME/MITgcm/build_options_shearwater" # build options file

# make directories
[ ! -d $build_dir ] && mkdir $build_dir
cd $build_dir

$mitgcm_root/tools/genmake2 -mpi -enable=mnc -mods "$mod_dir" -optfile "$buildopts" --rootdir "$mitgcm_root"
make -j 128 clean
make -j 128 depend
make -j 128
mv mitgcmuv $exec_dir/mitgcm_kwics.ex
