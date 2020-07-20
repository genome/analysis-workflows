#!/bin/bash

dir='/gscmnt/gc2698/jin810/references/split_regions'

cd $dir

interval_files=$( ls $PWD/*)

for i in ${interval_files[@]}
do
    echo - $i
done
