#!/usr/bin/env python
# -*- coding: utf-8 -*-]
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from math import log


"""
arguments:
	seqDict - sequence dictionary
"""


def consensus(seqDict):
	aaDict={}	#dictionary od aa's frequency
	id=seqDict.keys()	#seq IDs
	h=0		#uncertainty
	r=0		#information content
	height={}

	for i in xrange(len(seqDict[id[0]])):
		for key in id:
			if aaDict.has_key(seqDict[key][i]):
				aaDict[seqDict[key][i]]+=1.0
			else:
				aaDict[seqDict[key][i]]=1.0

		for k in sorted(aaDict):
			aaDict[k]/=len(seqDict)
			h+=aaDict[k]*log(aaDict[k], 2) # + zamiast - -, by nie robić niepotrzebnych błędów numeryccznych

		r=log(20,2)+h
		h=0

		for k in sorted(aaDict):
			height[k]=aaDict[k]*r
			print "position", i, "height", k, height[k]

		aaDict={}

#TODO: krytetia wyboru znakow do konsensusu