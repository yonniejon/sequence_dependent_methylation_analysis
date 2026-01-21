#!/bin/bash
zcat "$1" | awk -F"\t" -v OFS="\t" '{t=0.3; a=$5+$6+$7+1; b=$8+$9+$10+1; prop_a=a/(a+b); prop_b=b/(a+b); if((prop_a>t) && (prop_b>t) && (a > 5 && b > 5)){print $1, $2, "1"}else if(a > 10 && b < 2){print $1, $2, "2"}else if(a < 2 && b > 10){print $1, $2, "3"}else{print $1, $2, "0"}}' > "$2";
