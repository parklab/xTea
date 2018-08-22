import os

cmd="minimap2 -Hk19 -Xw3 -m40 -g10000 --max-chain-skip 25 -t8 21_29069186.fa 21_29069186.fa > 21_29069186.paf"
cmd="miniasm -s 100 -c 1 -e 1  -f 21_29069186.fa 21_29069186.paf > 21_29069186.untig"
cmd="minimap2 -Hk19 -Xw3 -m40 -g10000 --max-chain-skip 25 21_29069186_utg.fa  21_29069186.fa >  21_29069186_2_utg.paf"
cmd="racon 21_29069186.fa 21_29069186_2_utg.paf 21_29069186_utg.fa"

