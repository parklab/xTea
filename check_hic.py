import sys
import os 

sf_pos=sys.argv[1]

tei_pos=147020933
tei_start=147010933
tei_end=147030933
tei_chrm="chr2"
with open(sf_pos) as fin_pos:
	for line in fin_pos:
		fields=line.split()
		chrm=fields[0]
		pos=int(fields[1])
		pos_start=pos-10000
		pos_end=pos+10000

		ss="{0}:{1}-{2} & {3}:{4}-{5} [offset 0,0:0,0]".format(tei_chrm, tei_start, tei_end, chrm, pos_start, pos_end)
		print ss

