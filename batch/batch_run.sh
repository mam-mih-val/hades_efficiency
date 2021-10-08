#!/bin/bash

format='+%Y/%m/%d-%H:%M:%S'

date $format

job_num=$(($SLURM_ARRAY_TASK_ID))

filelist=$lists_dir/$(ls $lists_dir | sed "${job_num}q;d")

cd $output_dir
mkdir -p $job_num
cd $job_num

while read line; do
    echo $line >> list.txt
done < $filelist
echo >> list.txt

module use /cvmfs/vae.gsi.de/centos7/modules/linux-centos7-x86_64/cmake-3.17.3-gcc-8.1.0-3354fvg

echo "loading /lustre/hades/user/mmamaev/install/root-6.18.04-centos7-cxx17/bin/thisroot.sh"
source /lustre/hades/user/mmamaev/install/root-6.18.04-centos7-cxx17/bin/thisroot.sh

echo "executing $build_dir/acceptance list.txt $pdg_code"
$build_dir/acceptance list.txt $pdg_code

echo JOB FINISHED!