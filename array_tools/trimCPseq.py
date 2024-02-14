#!/usr/bin/env python
""" Trim an overly-long CPseq file
#Lena- I think this is for a chip that has a too high-density of RNA, but I'm not toally sure- asking...

 Inputs:
   Sequence data (.CPseq files)

 Outputs:
   trimmed Sequence data

 Ben Ober-Reynolds, boberrey@stanford.edu
 20160805

 """

import sys
import os
import argparse


### MAIN ###


def main():
	################ Parse input parameters ################

	#set up command line argument parser
	parser = argparse.ArgumentParser(description='Script for trimming CPseq files')
	group = parser.add_argument_group('required arguments for processing data')
	group.add_argument('-f', '--file', required=True,
	                    help='file to be trimmed (CPseq)')
	group.add_argument('-n', '--nth_line', required=True,
	                    help='Every nth line will be trimmed')

	group = parser.add_argument_group('optional arguments')
	# Currently this option doesn't do anything. May be worth implementing in the future
	#group.add_argument('-rf','--respect_filter', default=True,
	#                    help='do not trim clusters that have a filter')


	if not len(sys.argv) > 1:
	    parser.print_help()
	    sys.exit()

	args = parser.parse_args()

	seqFileName = args.file
	n = int(args.nth_line)
	if not os.path.isfile(seqFileName):
		print "File "+seqFileName+" is not a valid file! Exiting..."
		sys.exit()

	with open(seqFileName, 'r') as f:
		seqData = f.readlines()

	print "Trimming every "+str(n)+"th line from file "+seqFileName+"..."

	newFileName = "trimmed_"+str(n)+"_"+seqFileName
	with open(newFileName, 'w') as f:
		count = 1
		for line in seqData:
			fields = len(line.split())
			if count == n: 
				if fields > 9:
					# Don't throw away clusters with an assigned filter
					f.write(line)
				count = 1
			else:
				f.write(line)
				count += 1
	print "...done"








if __name__ == '__main__':
    main()
