#!/usr/bin/env python3

import sys
import re

def usage(message = None):
	if message is not None:
		print(message)
	print('Usage: readtrimmer <infile> <number of nucleotides 3 prime> <quality of nucleotides 3 prime> <number of nucleotides 5 prime> <quality of nucleotides 5 prime>')
	sys.exit(1)

if len(sys.argv) < 2:
	usage()


try:
	infile = open(sys.argv[1], 'r')
	outfile = open(str(sys.argv[1]) + '_trimmed', 'w')
except IOError:
	usage("Can't open file, " + sys.argv[1])


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

def detect_phred(file):
	phred33_unique = {'!', '"', '#', '$', '%', '&', "'", '(', ')', '*', 
	'+', ',', '-', '.', '/', '0', '1', '2', '3', '4',
	'5', '6', '7', '8', '9', ':', ';', '<', '=', '>',
	'?'}
	phred64_unique = {'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
	'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^',
	'_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
	'i'}
	phred_scale = ''
	flag = False
	for line in infile:
		if line.startswith('@HWI'):
			flag = False
		if flag is True:
			quality_data = line
			for asci in quality_data:
				if asci in phred33_unique:
					phred_scale = 33
					return phred_scale 
				elif asci in phred64_unique:
					phred_scale = 64
					return phred_scale
				else:
					break
		if line.startswith('+'):
			flag = True
		if phred_scale != '':
			break

print(detect_phred(infile))


seq = 'NGTACCAGCTGCTGCAAAATCATGCCAAAATGTAAAGACCATCGAGACTAGGAAGAAACTGCATCAACTAACGAGCAAAATAACCAGCTAACATCATAATG'
qual = '#11=ABDDHFDFHIIIIIIEHIC?@GHIIICG?CFHEHGIICF;BDEG9B;FGHHIIII4=C@=CEGH:CHFEE<99CCC@CCCCC@B5:ACC@AC>AC##'
seq_trimmed_5 = ''
seq_trimmed_3 = ''
asci_5 = ''
asci_3 = ''
window_size = 5



seq_qual_5 = list(zip(seq,qual))
print(seq_qual_5)
"""
#Trim from 5' based on minimum quality of single nucleotides
for nuc,asci in seq_qual_5:
	if phred33[asci]<30 and seq_trimmed_5 == '':
		continue
	else:
		seq_trimmed_5 += nuc
		asci_5 += asci

seq_qual_3 = list(zip(seq_trimmed_5,asci_5))

#Trim from 3' based on minimum quality of single nucleotides
for nuc,asci in seq_qual_3[::-1]:
	if phred33[asci]<30 and seq_trimmed_3 == '':
		continue
	else:
		seq_trimmed_3 += nuc
		asci_3 += asci
"""

#Trim from 5' based on average of moving window
i = 0
while i < len(seq_qual_5) - window_size + 1:
	sum_window = 0
	window = seq_qual_5[i:i + window_size]
	for nuc,asci in window:
		sum_window += phred33[asci]
	window_average = sum_window/window_size
	if window_average < 30 and seq_trimmed_5 == '':
		i += 1
		continue
	else:
		i += 1
		seq_trimmed_5 += window[0][0]
		asci_5 += window[0][1]

#Add remaining 3' nucleotides to trimmed sequence
for nuc,asci in window[1:]:
	seq_trimmed_5 += nuc
	asci_5 += asci

seq_qual_3 = list(zip(seq_trimmed_5,asci_5))
seq_qual_3 = seq_qual_3[::-1]
print(len(seq_qual_3))

#Trim from 3' based on average of moving window
i = 0
while i < len(seq_qual_3) - window_size + 1:
	sum_window = 0
	window = seq_qual_3[i:i + window_size]
	for nuc,asci in window:
		sum_window += phred33[asci]
	window_average = sum_window/window_size
	if window_average < 30 and seq_trimmed_3 == '':
		i += 1
		continue
	else:
		i += 1
		seq_trimmed_3 += window[0][0]
		asci_3 += window[0][1]

#Add remaining 5' nucleotides to trimmed sequence
for nuc,asci in window[1:]:
	seq_trimmed_3 += nuc
	asci_3 += asci

print(seq)
print(qual)
print(seq_trimmed_5)
print(asci_5)
print(seq_trimmed_3[::-1])
print(asci_3[::-1])


