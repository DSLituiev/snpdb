#!/bin/sh
python3.4 ../vcf_to_sqlite /media/Processing/seq/data/MO7-rmdup-clipOverlap-q20-freebayes-ems-annotation.vcf -l MO7

python3.4 ../vcf_to_sqlite /media/Processing/seq/data/MO8-rmdup-clipOverlap-q20-freebayes-ems-annotation.vcf -l MO8
