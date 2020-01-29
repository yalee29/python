import pysam

if __name__ == '__main__':

	import argparse

	parser = argparse.ArgumentParser()
	parser.add_argument('bam')
	args = parser.parse_args()

	bam = pysam.AlignmentFile(args.bam, "rb")
	
	line = 0
	
	for read in bam.fetch('MT',0,16569 ):

		cigar = read.cigartuples
		if cigar is  None:
			continue
		
		if cigar[0][0] ==4 or cigar[len(cigar)-1][0] == 4:
			try:
				print(read.query_name, read.reference_start, read.reference_end)
				line += 1
			except:
				continue
#	print(line)
