intersectBed -a $1 -b $2 -f 0.5 -wo > overlap.temp
intersectBed -a $2 -b $1 -f 0.5 -wo >> overlap.temp
intersectBed -a $1 -b $3 -f 0.5 -wo >> overlap.temp
intersectBed -a $3 -b $1 -f 0.5 -wo >> overlap.temp
intersectBed -a $2 -b $3 -f 0.5 -wo >> overlap.temp
intersectBed -a $3 -b $2 -f 0.5 -wo >> overlap.temp
cut -f 1-3 overlap.temp > to_merge.temp
cut -f 11-13 overlap.temp >> to_merge.temp
sort -k1,1 -k2,2n to_merge.temp | mergeBed -i - > $4
