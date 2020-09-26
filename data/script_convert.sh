#!/bin/bash

declare -a arr=(
    "Run3v2_Barrel_NThltIter2_0"
    "Run3v2_Barrel_NThltIter2_1"
    "Run3v2_Barrel_NThltIter2_2"
    "Run3v2_Barrel_NThltIter2_3"
    "Run3v2_Barrel_NThltIter2FromL1_0"
    "Run3v2_Barrel_NThltIter2FromL1_1"
    "Run3v2_Barrel_NThltIter2FromL1_2"
    "Run3v2_Barrel_NThltIter2FromL1_3"
    "Run3v2_Endcap_NThltIter2_0"
    "Run3v2_Endcap_NThltIter2_1"
    "Run3v2_Endcap_NThltIter2_2"
    "Run3v2_Endcap_NThltIter2_3"
    "Run3v2_Endcap_NThltIter2FromL1_0"
    "Run3v2_Endcap_NThltIter2FromL1_1"
    "Run3v2_Endcap_NThltIter2FromL1_2"
    "Run3v2_Endcap_NThltIter2FromL1_3"
)

for i in "${arr[@]}"
do
  strin=$i
  strout=${strin/NThlt/hlt}
  xmllint --format "$strin.xml" > "$strout.xml"
  rm "$strin.xml"
done;

echo "finish"
