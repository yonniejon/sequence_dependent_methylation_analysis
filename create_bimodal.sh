#!/bin/bash
zcat "$1" | awk 'function min(a,b) {return a < b ? a : b} { weakest=min($6, $8);total_c=$6 + $7 + $8; min_u_m=total_c > 0 ? weakest/($6 + $7 + $8) : 0; if((min_u_m > 0.2) && (total_c > 25)){print $0}}' > "$2".bimodal.tsv
