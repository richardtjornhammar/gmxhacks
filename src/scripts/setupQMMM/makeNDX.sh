#!/bin/bash
touch made.ndx
echo "[ QMatoms ]"> made.ndx
for k in `seq 1 1 30`;
do 
echo $k >> made.ndx
done

 

