#!/usr/bin/env python3
import sys
import re

#### 2. Detect Phred score ####
infile=open("sequence_example.txt","r")


def translation(quality_data):
    '''translates a line of quality scores to a numerical value'''
    quality_table = {
	'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8, '*': 9, 
	'+': 10, ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19,
	'5': 20, '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29,
	'?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39,
	'I': 40, 'J': 41, 'K': 42, 'L': 43, 'M': 44, 'N': 45, 'O': 46, 'P': 47, 'Q': 48, 'R': 49,
    'S': 50, 'T': 51, 'U': 52, 'V': 53, 'W': 54, 'X': 55, 'Y': 56, 'Z': 57, '[': 58, '\\': 59,
    ']': 60, '^': 61, '_': 62, '`': 63, 'a': 64, 'b': 65, 'c': 66, 'd': 67, 'e': 68, 'f': 69,
    'g': 70, 'h': 71}
    #make data to a list
    data_list = []
    for char in quality_data:
        data_list.append(char)
    quality_scores = []
    for char in range(len(data_list)):
        quality_scores += [quality_table[char]]
    return(quality_scores)


def detect_phred(infile):
    '''detects phred score for sequences in file'''
    quality_data = None
    line_count = 0
    #iterate over file until scale is identified
    while phred_scale == None:
        for line in infile:
            line_count += 1    
            #identify quality score lines (every 4th)
            if line_count % 4 == 0:
                #translate the line
                line = translation(line)
                #identify if in phred 33 range
                print(line) 


#Set start
phred_scale = None 
         
print(detect_phred(infile))





