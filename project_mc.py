#!/usr/bin/env python3

import sys
import re
import argparse


####argument parser####
def check_pos(value):
    '''validate if args is a pos int'''
    try:
        int_value = int(value)
        if int_value < 0:
            raise argparse.ArgumentError("%s is not a positive int value" % value)
    except argparse.ArgumentError as error:
        print("%s is not a positive int value" % value)
        sys.exit()
    return int_value

def check_qual_cut(value):
    '''check if value is pos int or return defaultvalue if yes'''
    if value == 'yes':
        value = 0
    else:
        try:
            int_value = int(value)
            if int_value < 0:
                raise argparse.ArgumentError("%s is not a positive int value" % value)
        except argparse.ArgumentError as error:
            print("%s is not a positive int value" % value)
            sys.exit()
    return int_value

parser = argparse.ArgumentParser(description = 'Trimming if fastq file.')
parser.add_argument('-outfilename', required = True, 
                    help = 'output filename. If gzip file is wanted, add .gz in end of filename')
parser.add_argument('-filename', required = True, 
                    help = 'input fastq filename for trimming')
parser.add_argument('-phred_scale', choices = ['phred+33', 'phred+64'], 
                    help = 'phred scale (type: phred+33 or phred+64') 
parser.add_argument('-fixed_trim', nargs = 2, default=[0,0], type = check_pos, 
                    help = 'what fixed base length to trim from each end (type: space seperated string of pos int of length 2')
parser.add_argument('-min_residue3', type = check_qual_cut, default = False, 
                    help = 'if trim from 3end should be minimum of single residue and what quality should be cutoff value (type: pos int / yes, default cutoff = 40)') 
parser.add_argument('-min_residue5', type = check_qual_cut, default = False, 
                    help = 'if trim from 5end should be minimum of single residue and what quality should be cutoff value (type: pos int / yes, default cutoff = 40)') 
parser.add_argument('-mean_mw3', type = check_qual_cut, default = False, 
                    help = 'if trim from 3end should be mean of moving window and what cutoff mean should be (type: pos int / yes, default: 40)')   
parser.add_argument('-mean_mw5', type = check_qual_cut, default = False, 
                    help = 'if trim from 5end should be mean of moving window and what cutoff mean should be (type: pos int / yes, default: 40)') 

args = parser.parse_args()

#raise exception if moving window and single residue trimming methods is used for same ends
if args.min_residue3 != False and args.mean_mw3 != False:
    raise argparse.ArgumentError('Cannot trim with both moving window and single residue from 3end')
if args.min_residue5 != False and args.mean_mw5 != False:
    raise argparse.ArgumentError('Cannot trim with both moving window and single residue from 5end')    

#set default trimming if nothing is specified



####1. read file####

#read input file if filetype is gzipped fastq 
if args.filename.endswith('fastq.gz'):
    try:
        infile = gzip.open(args.filename,'rt')
    except IOError as error:
        sys.stdout.write('Cannot open inputfile, reason: ' + str(error) + '\n')
        sys.exit(1)
#read file if filetype is fastq
elif args.filename.endswith('fastq'):
    try:
        infile = open(args.filename, 'r')
    except IOError as error:
        sys.stdout.write('Cannot open inputfile, reason: ' + str(error) + '\n')
        sys.exit(1)
else:
    print('Not a valid filename')

#open output file and check if name if gzip
if args.outfilename.endswith('gz'):
    try:
        outfile = gzip.open(args.outfilename, 'wt')
    except IOError as error:
        sys.stdout.write('Cannot open outfile, reason: ' + str(error) + '\n')
        sys.exit(1)
else:
    try:
        outfile = open(args.outfilename, 'w')
    except IOError as error:
        sys.stdout.write('Cannot open outfile, reason: ' + str(error) + '\n')
        sys.exit(1)   




def usage(message = None):
	if message is not None:
		print(message)
	print('Usage: readtrimmer <infile> <number of nucleotides 3 prime> <quality of nucleotides 3 prime> <number of nucleotides 5 prime> <quality of nucleotides 5 prime>')
	sys.exit(1)


#### 2. Detect Phred score ####

phred33 = {
	'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8, '*': 9, 
	'+': 10, ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19,
	'5': 20, '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29,
	'?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39,
	'I': 40, 'J': 41}
phred64 = {
	'@': 0, 'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8,	'I': 9,
	'J': 10, 'K': 11, 'L': 12, 'M': 13, 'N': 14, 'O': 15, 'P': 16, 'Q': 17, 'R': 18, 'S': 19,
	'T': 20, 'U': 21, 'V': 22, 'W': 23, 'X': 24, 'Y': 35, 'Z': 26, '[': 27, '\\': 28, ']': 29,
	'^': 30, '_': 31, '`': 32, 'a': 33, 'b': 34, 'c': 35, 'd': 36, 'e': 37, 'f': 38, 'g': 39, 
	'h': 40, 'i': 41}



def detect_phred(infile):
    '''detects phred score for sequences in file'''
    #Set start
    phred_scale = None 
    quality_data = None
    line_count = 0
    #iterate over file until scale is identified
    while phred_scale == None:
        for line in infile:
            line = line.strip()
            line_count += 1    
            #identify quality score lines (every 4th)
            if line_count % 4 == 0:
                for char in line:
                    #translate the line to find phred scales
                    if char in phred33 and char not in phred64:
                        phred_scale = "phred+33"
                    elif char in phred64 and char not in phred33:
                        phred_scale = "phred+64"
    #if program can not determine phred scale
    if phred_scale == None:
        print("Cannot determine phred scale")
    return phred_scale    
                  

def trim_fixed(seq_qual, nuc_to_trim_5, nuc_to_trim_3):
	'''Trim fixed number of nucleotides from 5' and 3' ends'''
	count = 0
	seq_trimmed_3 = ''
	asci_trimmed_3 = ''
	seq_trimmed_5 = ''
	asci_trimmed_5 = ''
	for nuc, asci in seq_qual:
		if count < nuc_to_trim_5:
			count += 1
		else:
			seq_trimmed_5 += nuc
			asci_trimmed_5 += asci
	seq_qual_trimmed_5 = list(zip(seq_trimmed_5, asci_trimmed_5))
	count = 0
	for nuc, asci in seq_qual_trimmed_5[::-1]:
		if count < nuc_to_trim_3:
			count += 1
		else:
			seq_trimmed_3 += nuc
			asci_trimmed_3 += asci
	seq_qual_trimmed_3 = list(zip(seq_trimmed_3[::-1], asci_trimmed_3[::-1]))
	return seq_qual_trimmed_3


def trim_single_nuc_5(seq_qual, threshold):
	'''Trim from 5' based on minimum quality of single nucleotides'''
	seq_trimmed_5 = ''
	asci_trimmed_5 = ''
	for nuc,asci in seq_qual:
		if phred_scale == 'phred+33':
			if phred33[asci] < threshold and seq_trimmed_5 == '':
				continue
			else:
				seq_trimmed_5 += nuc
				asci_trimmed_5 += asci
		elif phred_scale == 'phred+64':
			if phred64[asci] < threshold and seq_trimmed_5 == '':
				continue
			else:
				seq_trimmed_5 += nuc
				asci_trimmed_5 += asci
		else:
			raise ValueError('Phred scale not correctly identified', phred_scale)
			sys.exit(1)
	seq_qual = list(zip(seq_trimmed_5, asci_trimmed_5))
	return seq_qual


def trim_single_nuc_3(seq_qual, threshold):
	'''Trim from 3' based on minimum quality of single nucleotides'''
	seq_trimmed_3 = ''
	asci_trimmed_3 = ''
	for nuc,asci in seq_qual[::-1]:
		if phred_scale == 'phred+33':
			if phred33[asci] < threshold and seq_trimmed_3 == '':
				continue
			else:
				seq_trimmed_3 += nuc
				asci_trimmed_3 += asci
		elif phred_scale == 'phred+64':			
			if phred64[asci] < threshold and seq_trimmed_3 == '':
				continue
			else:
				seq_trimmed_3 += nuc
				asci_trimmed_3 += asci	
		else:
			raise ValueError('Phred scale not correctly identified', phred_scale)
			sys.exit(1)
	seq_qual = seq_trimmed_3[::-1], asci_trimmed_3[::-1]
	return seq_qual


def trim_moving_window_5(seq_qual, threshold):
	'''Trim from 5' based on average of moving window of size 5'''
	seq_trimmed_5 = ''
	asci_trimmed_5 = ''
	window_size = 5
	i = 0
	while i < len(seq_qual) - window_size + 1:
		sum_window = 0
		window = seq_qual[i:i + window_size]
		if phred_scale == 'phred+33':
			for nuc,asci in window:
				sum_window += phred33[asci]
			window_average = sum_window/window_size
			if window_average < threshold and seq_trimmed_5 == '':
				i += 1
				continue
			else:
				i += 1
				seq_trimmed_5 += window[0][0]
				asci_trimmed_5 += window[0][1]
		elif phred_scale == 'phred+64':
			for nuc,asci in window:
				sum_window += phred64[asci]
			window_average = sum_window/window_size
			if window_average < threshold and seq_trimmed_5 == '':
				i += 1
				continue
			else:
				i += 1
				seq_trimmed_5 += window[0][0]
				asci_trimmed_5 += window[0][1]
		else:
			raise ValueError('Phred scale not correctly identified', phred_scale)
			sys.exit(1)

	#Add remaining 3' nucleotides to trimmed sequence
	for nuc,asci in window[1:]:
		seq_trimmed_5 += nuc
		asci_trimmed_5 += asci

	seq_qual = list(zip(seq_trimmed_5, asci_trimmed_5))
	return seq_qual


def trim_moving_window_3(seq_qual, threshold):
	'''Trim from 3' based on average of moving window of size 5'''
	seq_trimmed_3 = ''
	asci_trimmed_3 = ''
	window_size = 5
	i = 0
	seq_qual = seq_qual[::-1]
	while i < len(seq_qual) - window_size + 1:
		sum_window = 0
		window = seq_qual[i:i + window_size]
		if phred_scale == 'phred+33':
			for nuc,asci in window:
				sum_window += phred33[asci]
			window_average = sum_window/window_size
			if window_average < threshold and seq_trimmed_3 == '':
				i += 1
				continue
			else:
				i += 1
				seq_trimmed_3 += window[0][0]
				asci_trimmed_3 += window[0][1]
		elif phred_scale == 'phred+64':
			for nuc,asci in window:
				sum_window += phred64[asci]
			window_average = sum_window/window_size
			if window_average < threshold and seq_trimmed_3 == '':
				i += 1
				continue
			else:
				i += 1
				seq_trimmed_3 += window[0][0]
				asci_trimmed_3 += window[0][1]
		else:
			raise ValueError('Phred scale not correctly identified', phred_scale)
			sys.exit(1)

	#Add remaining 5' nucleotides to trimmed sequence
	for nuc,asci in window[1:]:
		seq_trimmed_3 += nuc
		asci_trimmed_3 += asci

	seq_qual = seq_trimmed_3[::-1], asci_trimmed_3[::-1]
	return seq_qual


#MAIN PROGRAM

#Initialize
read_count = 0
line_count = 0
seq = ''
qual = ''


phred_scale = detect_phred(infile)

infile.seek(0)
for line in infile:
	line = line.strip()
	line_count += 1
	#Identify header, sequence and quality data for read
	if line_count % 4 == 1:
		header = line
	elif line_count % 4 == 2:
		seq = line
	elif line_count % 4 == 0:
		qual = line
	#If sequence and quality data are both containing data - trimming functions are performed according to arguments from command line
	if seq != '' and qual != '':
		seq_qual = list(zip(seq, qual))
		#Trim fixed number of nucleotides from 5' end
		seq_qual_trimmed_fixed = trim_fixed(seq_qual, args.fixed_trim[0], args.fixed_trim[0])
		#Trim based on quality of nucleotides from 5' end
		if args.min_residue5 is not False:
			seq_qual_trim5 = trim_single_nuc_5(seq_qual_trimmed_fixed, args.min_residue5)
		else:
			seq_qual_trim5 = trim_moving_window_5(seq_qual_trimmed_fixed, args.mean_mw5)
		#Trim based on quality of nucleotides from 3' end
		if args.min_residue3 is not False:
			seq_qual_trim3 = trim_single_nuc_3(seq_qual_trim5, args.min_residue3)
		else:
			seq_qual_trim3 = trim_moving_window_3(seq_qual_trim5, args.mean_mw3)

		#Print read to outfile
		print(header, file = outfile)
		print(seq_qual_trim3[0], file = outfile)
		print('+', file = outfile)
		print(seq_qual_trim3[1], file = outfile)

		#Increase read count by 1 and reset seq and qual for next read
		read_count += 1
		seq = ''
		qual = ''

print('Number of reads in', args.filename, ':', read_count)
