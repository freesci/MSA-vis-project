#!/usr/bin/env python
# -*- coding: utf-8 -*-]

import textwrap
import os.path
import sys

try:
	import argparse
except ImportError:
	print "Argparse packages not found. You may download it from: http://argparse.googlecode.com/files/argparse-1.2.1.tar.gz"
	pass

try:
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio import SeqIO
except ImportError:
	print "Biopython packages not found. You may download it from: http://www.biosino.org/mirror/www.biopython.org/Download/default.htm"
	pass

try:
	from termcolor import colored, cprint
except ImportError:
	print "Termcolor packages not found. You may download it from: http://pypi.python.org/pypi/termcolor"
	pass

def inpout():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Visualisation of MSA',prog="MSA",
                                     epilog=textwrap.dedent('''\
                                     Details:
                                     ----------
                                     File with the MSA should contain sequences of the same length,
                                     preceded by an identifier, which occurs after the ">".
                                     Repetitive sequences will be omitted.
                                     In other cases the file will be regarded as inpcorrect.'''))

    parser.add_argument('inp', type=str, nargs='*', help='file .fasta with MSA')

    parser.add_argument('-o', '--output', type=str, nargs='?',default="MSAvis",
                        help='the name of the file where MSA visualization will be saved')

    parser.add_argument('-a', '--aminoacids', type=int, nargs='?',default=30,
                        help='number of aminoacids in one row in graph. Enter the number greater than 0 (Default 30).')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0\nAuthors: M.Habich, M.Maksymiuk, M.Stepniewska, K.Wreczycka')

    args = parser.parse_args()
    if len(os.path.splitext(args.output))>1:
        args.output=os.path.splitext(args.output)[0]

    if args.aminoacids < 1:
        args.aminoacids=30

    return args.inp, args.output, args.aminoacids, parser

inp,out,nam,parser=inpout()

def checkinp(parser,inp):
    if len(inp)==0:
        print colored("Error: Not specified data file\n----------------------------------",'red')
        parser.print_help()
        sys.exit()
    return inp[0]

inp=checkinp(parser, inp)

def seqDict(ifile):
    def has_key(d,key):
        if d.has_key(key.id):
            if d[key.id]!=key.seq.data:
                key.id=key.id+"(1)"
    	return key

    def equal_length(seq,num,lens):
        if len(seq.seq.data)!=lens:
            print colored("Error: Sequences number "+str(num)+" aren't of the same length like the other\n-----------------------------------------",'red')
            parser.print_help()
            sys.exit()

    seqDict={}
    iteration=1

    try:
	handle = open(ifile, "rU")
    except IOError:
	print colored("Error: File '"+ifile+"' doesn't exist\n----------------------------",'red')
        sys.exit()

    try:
        for record in SeqIO.parse(handle, "fasta") :
            if len(seqDict)==0: 
                seqDict[record.id]=record.seq.data
	        lenseq=len(record.seq.data)
            else:
                record=has_key(seqDict,record)
		equal_length(record,iteration,lenseq)
                seqDict[record.id]=record.seq.data
		lenseq=len(record.seq.data)
	    iteration+=1
    except IndexError:
        print colored("Error: Sequence without ID found\n---------------------------------------------",'red')
        parser.print_help()
        sys.exit()

    handle.close()

    if len(seqDict)==0:
        print colored("Error: The specified file is empty\n----------------------------",'red')
        parser.print_help()
        sys.exit()     

    return seqDict

dictionary=seqDict(inp)
if nam>len(dictionary):
	nam=len(dictionary)

