#!/bin/bash

files=*.pdf

for file in $files
do
  pdf2ps $file
done  


files=*.ps
for file in $files
do
  ps2eps $file
done

rm *.ps  
rm *.pdf

#for i in `ls *.eps`
#do
#  convert $i ${i%.*}.jpg
#done

#rm *.eps

