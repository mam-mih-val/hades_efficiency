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

echo "executing $build_dir/acceptance list.txt $pdg_code /lustre/hebe/hades/user/bkardan/param/centrality_epcorr_mar19_ag123ag_2500A_glauber_gen4_2020_10_pass3.root"
$build_dir/acceptance list.txt $pdg_code /lustre/hebe/hades/user/bkardan/param/centrality_epcorr_mar19_ag123ag_2500A_glauber_gen4_2020_10_pass3.root

echo JOB FINISHED!