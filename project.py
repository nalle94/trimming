#!/usr/bin/env python3
import sys
import re
import gzip
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
    '''check if value is pos int or return defaultvalue if not'''
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

parser = argparse.ArgumentParser(description = 'Trimming of fastq file.')
parser.add_argument('-o', dest='outfilename', required = True, 
                    help = 'output filename. If gzip file is wanted, add .gz in end of filename')
parser.add_argument('-f', dest='filename', required = True, 
                    help = 'input fastq filename for trimming')
parser.add_argument('-p', dest='phred_scale', choices = ['phred+33', 'phred+64'], 
                    help = 'phred scale (type: phred+33 or phred+64)') 
parser.add_argument('-t', dest='fixed_trim', nargs = 2, default=[0,0], type = check_pos, 
                    help = 'what fixed base length to trim from each end (type: space seperated string of pos int of length 2')
parser.add_argument('-m3', dest='min_residue3', type = check_qual_cut, default = False, 
                    help = 'if trim from 3end should be minimum of sigle residue and what quality should be cutoff value (type: pos int / yes, default cutoff = 40)') 
parser.add_argument('-m5', dest='min_residue5', type = check_qual_cut, default = False, 
                    help = 'if trim from 5end should be minimum of sigle residue and what quality should be cutoff value (type: pos int / yes, default cutoff = 40)') 
parser.add_argument('-w3', dest='mean_mw3', type = check_qual_cut, default = False, 
                    help = 'if trim from 3end should be mean of moving window and what cutoff mean should be (type: pos int / yes, default: 40)')   
parser.add_argument('-w5', dest='mean_mw5', type = check_qual_cut, default = False, 
                    help = 'if trim from 5end should be mean of moving window and what cutoff mean should be (type: pos int / yes, default: 40)') 

args = parser.parse_args()

#raise exception if moving window and single residue trimming methods is used for same ends
if parser.min_residue3 != False or parser.mean_mw3 != False:
    raise argparse.ArgumentError('Cannot trim with both moving window and single residue from 3end')
if parser.min_residue5 != False or parser.mean_mw5 != False:
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


#### 2. Detect Phred score ####
infile=open("sequence_example.txt","rt")
outfile = gzip.open("sequence_example.txt.gzip","wt")

infile.close()
outfile.close()


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
            





