#!/bin/bash
head -n 1 $1

for f in $*
do
     tail -n +2 ${f}
done
