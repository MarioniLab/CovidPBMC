#! /usr/bin/bash

## nmove files and rename them according to the input mapping

while IFS= read -r line;
do
    #echo $line
    ln -s $line
done < $1
