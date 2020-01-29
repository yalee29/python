import pysam
import pandas as pd
import argparse
import re

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('bam')
	parser.add_argument('-o', '--out', help="output tsv file")
	parser.add_argument('-d', '--debug', help="debug mode", action='store_true')
	args = parser.parse_args()

	bam = pysam.AlignmentFile(args.bam, "rb")
	
	d = {}
	for read in bam.fetch('MT',0,16569 ):

		cigar = read.cigartuples

		if cigar is  None:
			continue
		
		try:
			SA = read.get_tag("SA")
		except:
			continue
		

		sup_chr = SA.split(',')[0]
		sup_pos = SA.split(',')[1]
		sup_cigar = SA.split(',')[3]

		
		if sup_chr != "MT":
			continue
		

		if cigar[0][0] ==4 or cigar[len(cigar)-1][0] == 4:

			if cigar[0][0] ==4 and (cigar[len(cigar)-1][0] != 4):
				bp1 = read.reference_start
			elif cigar[0][0] !=4 and (cigar[len(cigar)-1][0] == 4):
				bp1 = read.reference_end
			else:
				continue  #일단 양쪽 말단에 soft clip이 있을 경우는 제외
		
			p1 = re.compile('^\d+S\d+M$')
			p2 = re.compile('^\d+M\d+S$')
			if p1.match(sup_cigar): 
				bp2 = sup_pos
			elif p2.match(sup_cigar):
				bp2 = int(sup_pos) + int(sup_cigar.split('M')[0]) - 1
			else:
				continue
		else:
			continue

		bp = str(bp1) + '_' + str(bp2)

		if args.debug:
			print("Read name:", read.query_name,"Break points:", bp, "Mapping quality:",read.mapping_quality )
	
		if bp in d.keys():
			d[bp]+=1
		else:
			d[bp]=1



	df =pd.DataFrame([{'break_point':b, 'n_reads':r} for (b,r) in d.items()])
	print(df)
	
	if not args.debug:
		df.to_csv(args.out, sep="\t", index=False)


