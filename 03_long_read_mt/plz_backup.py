import pysam

if __name__ == '__main__':

	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument('bam')
	parser.add_argument('-o', '--out', help="output tsv file", required=True)
	args = parser.parse_args()

	bam = pysam.AlignmentFile(args.bam, "rb")
	
#	line = 0
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

		
		if sup_chr == "MT":

			if cigar[0][0] ==4 or cigar[len(cigar)-1][0] == 4:
				print(read.query_name, read.reference_start, cigar, SA)
			#try:
				if cigar[0][0] ==4 and cigar[len(cigar)-1][0] == 4:
					if "5_"+ str(read.reference_start) + "_3_" + str(read.reference_end) in d.keys():
						d["5_"+ str(read.reference_start)+ "_3_"+ str(read.reference_end)] += 1
					else:
						d["5_"+ str(read.reference_start) + "_3_"+ str(read.reference_end)] = 1
#					print(read.query_name,"5'&3'", read.reference_start), read.reference_end
				elif cigar[0][0] ==4 and  (cigar[len(cigar)-1][0] != 4):
					if "5_"+ str(read.reference_start) in d.keys():
						d["5_"+ str(read.reference_start)] += 1
					else:
						d["5_"+ str(read.reference_start)] = 1
#					print(read.query_name,"5",read.reference_start))
				else:
					if "3_"+ str(read.reference_end) in d.keys():
						d["3_"+ str(read.reference_end)] += 1
					else:
						d["3_"+ str(read.reference_end)] = 1 
#					print(read.query_name,"3'" ,read.reference_end
		else:
			continue
			#	line += 1
			#except:
				#continue
	#	sup_chr = SA.split('.')[0]
	#	sup_pos = SA.split(',')[1]
	#	sup_cigar = SA.split(',')[3]
	#	print(sup_pos)
	#	print(sup_cigar)
#	print(line)

	import pandas as pd
	df =pd.DataFrame([{'break_point':b, 'n_reads':r} for (b,r) in d.items()])
#	print(df)
#	df.to_csv(args.out, sep="\t", index)


