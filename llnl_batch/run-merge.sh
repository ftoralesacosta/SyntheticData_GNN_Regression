#! /bin/bash
date
cd $1
echo 'Merging'
hadd total.root output/*root
date
