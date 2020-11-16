#!/usr/bin/env python3
import sys, re, gzip, argparse


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


parser = argparse.ArgumentParser(description = 'Trimming of fastq file.')
parser.add_argument('-f', dest = 'filename', required = True, 
                    help = 'input fastq filename for trimming')
parser.add_argument('-o', dest = 'outfilename', required = True, 
                    help = 'output filename. If gzip file is wanted, add .gz in end of filename')
parser.add_argument('-p', dest = 'phred_scale', choices = ['phred+33', 'phred+64'], 
                    help = 'phred scale (type: phred+33 or phred+64)') 
                    
#add grouped argument for trimming settings
trim = parser.add_argument_group('Trimming', description = 'Settings for trimming. Besides fixed_trim, only one trimming option can be performed from each end. If no options are entered, no trimming will be performed')
trim.add_argument('-t', dest = 'fixed_trim', nargs = 2, default = [0,0], type = check_pos, 
                    help = 'what fixed base length to trim from eachs end (type: space seperated string of pos int of length 2')
#add mutually exclusive argument for 3' end trimming
trim3 = trim.add_mutually_exclusive_group()
trim3.add_argument('-m3', dest = 'min_residue3', type = check_pos, default = False, 
                    help = 'minimum of single residue trimming from 3end (type: pos int, default: False)') 
trim3.add_argument('-a3', dest = 'mean_mw3', type = check_pos, default = False, 
                    help = 'mean of moving window trimming from 3end (type: pos int, default: False)')   
trim3.add_argument('-w3', dest = 'min_mw3', type = check_pos, default = False,
                    help = 'minimum of moving window trimming from 3end (type: pos int, default: False)')
#add mutually exclusive argument for 5' end trimming
trim5 = trim.add_mutually_exclusive_group()
trim5.add_argument('-m5', dest = 'min_residue5', type = check_pos, default = False, 
                    help = 'minimum of single residue trimming from 5end (type: pos int, default: False)') 
trim5.add_argument('-a5', dest = 'mean_mw5', type = check_pos, default = False, 
                    help = 'mean of moving window trimming from 5end (type: pos int, default: False)') 
trim3.add_argument('-w5', dest = 'min_mw5', type = check_pos, default = False,
                    help = 'minimum of moving window trimming from 5end (type: pos int, default: False)')

#add grouped argument for filtering settings
filtering = parser.add_argument_group('Filtering', description = 'Settings for filtering (default: min_mean_qual = 30, min_read_len = 50, max_N = 5)')
filtering.add_argument('-q', dest = 'min_mean_qual', type = check_pos, default = 30,
                    help = 'minimum mean quality after trimming (default = 30)')
filtering.add_argument('-r', dest = 'min_read_qual', type = check_pos, default = 50,
                    help = 'minimum read length after trimming (default = 50)')
filtering.add_argument('-s', dest = 'max_N', type = check_pos, default = 5,
                    help = 'maximum unknown N after trimming (default = 5)')

args = parser.parse_args()


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

#open log file
logfile = open('log.txt', 'w')


#### 2. Detect Phred score ####
infile=open("sequence_example.txt","rt")
outfile = gzip.open("sequence_example.txt.gzip","wt")



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
            

####close files####
infile.close()
outfile.close()
logfile.close()

