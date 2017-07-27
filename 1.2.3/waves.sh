#!/bin/bash

#pre-parsing (run only once, then comment)
python scripts/parser.py -i data-wos -o data -v

for y1 in {1970..2009} ;
 do mkdir -p data$y1;
 y2=$(($y1+3));
 echo $y1 $y2

 #step1: slice 4-year windows
 python  scripts/filter.py -i dat_files_clean -o data$y1 -y1 $y1 -y2 $y2 -v

 #step2: louvain method community detection (+gdf)
 python scripts/BC_mm_atomes.py -i data$y1 -o data$y1 -v

 # step3: create pdfs
 cd data$y1
 pdflatex *tex
 cd ..

done

