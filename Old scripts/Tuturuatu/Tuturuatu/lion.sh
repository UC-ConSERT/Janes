#!/bin/bash -e

awk '/frog/ {sum +=$3;nlines =NR} END {print sum/nlines}' test1.txt

pre_nlines=$(awk '/frog/ {print NR}' test1.txt)

wait

awk '/cat/ {sum +=$3;nlines =(NR-pre_nlines)} END {print sum/nlines}' test1.txt

wait

pre_nlines=$(awk '/cat/ {print NR}' test1.txt)

wait

awk '/lion/ {sum +=$3;nlines =(NR-pre_nlines)} END {print sum/nlines}' test1.txt

#######################
#!/bin/bash -e

frogsum=$(awk '/frog/ {sum +=$3} END {print sum}' test1.txt)

frogcount=$(grep -c frog test1.txt)

catmean=$(($(awk '/cat/ {sum +=$3} END {print sum}' test1.txt)/$(grep -c cat test1.txt)))  ##this works!!

catcount=$(grep -c cat test1.txt)

lionsum=$(awk '/lion/ {sum +=$3} END {print sum}' test1.txt)

lioncount=$(grep -c lion test1.txt)

lionmean=$((${frogsum}/${frogcount}))  ##this also works, but is longer!

echo "${lionmean}"
echo "${catmean}"