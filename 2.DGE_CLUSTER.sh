#!/bin/bash
lib=./lib/
dge_file=$1
myname=$2
$lib/plot.cluster.from.dge_file.r $dge_file $myname $cutoff
