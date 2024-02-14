#!/bin/bash
#
# Quantify array images on sherlock
#
# Usage: quantifyImages.sh prefix seq_dir roff_dir fluor_dir gv_path script_dir
#
# Ben Ober-Reynolds, boberrey@stanford.edu, 20191105

# prefix is the common prefix shared by each set of images (will likely be 'set')
prefix=$1
seq_dir=$2
roff_dir=$3
fluor_dir=$4
gv_path=$5
script_dir=$6

for s in *$prefix*
do
	echo $s
    mkdir $roff_dir/$s/
    mkdir $fluor_dir/$s/
    sbatch $script_dir/quantify_images.sbatch $s $seq_dir $roff_dir/$s/ $fluor_dir/$s/ $gv_path
done
