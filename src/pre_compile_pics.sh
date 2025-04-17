#! /bin/bash

for d in img-gen build; do
    if [ ! -d ${d} ]; then
    mkdir $d
    fi
done

for file in tikz/*.tex; do
    file_base=`basename $file`;

    cat img_IEEE.tex | sed -e "s/FILE/$file_base/" | \
    lualatex -output-directory build/ -jobname `basename $file_base .tex` && \
    mv build/`basename $file_base .tex`.pdf img-gen/
done
