#!/bin/bash
FILE=$1
while read LINE; do
     sh -c "$LINE"
done < $FILE
