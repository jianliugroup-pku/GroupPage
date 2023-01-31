#!/bin/bash

f=$1
fbase=${f%.*}
echo $f $fbase

pandoc $f -o $fbase.html --template=template.html
