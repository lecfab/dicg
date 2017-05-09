#!/bin/bash
for f in ../database/*.pgm ; do
    name=${f##*/}
    name=$(echo $name | sed 's/-.*pgm//g')
	echo "${name},$(./invariants $f)"
done