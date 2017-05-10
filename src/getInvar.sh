#!/bin/bash
echo "category,a,b,c,d,e,f,g,h,i"
for f in database/*.pgm ; do
    name=${f##*/}
    name=$(echo $name | sed 's/-.*pgm//g')
	echo "${name},$(./src/invariants $f)"
done