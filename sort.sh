#!/bin/bash

grep '@' $1 > sorted_temp.sam
grep -v '@' $1 | awk -F'[:\t]' '{print $8"\t"$10"\t"$11"\t"$0}' | sort -k1,1 -k2,2n -k3,3n | cut -f 4- >> sorted_temp.sam