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
	n=0		#number of non-gap letters at position i
	consensus=[]

	for i in xrange(len(seqDict[id[0]])):
		for key in id:
			if seqDict[key][i]!="-":
				if aaDict.has_key(seqDict[key][i]):
					aaDict[seqDict[key][i]]+=1.0
				else:
					aaDict[seqDict[key][i]]=1.0
				n+=1

		for k in sorted(aaDict):
			aaDict[k]/=n
			h+=aaDict[k]*log(aaDict[k], 2)

		r=log(20,2)+h
		h=0
		n=0

		for k in sorted(aaDict):
			height[k]=aaDict[k]*r


		first=("", 0)
		second=("", 0)
		for k in sorted(height):
			if height[k] > first[1]:
				first, second=(k, height[k]), first
			elif height[k] > second[1]:
				second=(k, height[k])

		if first[1]>log(20,2)*0.5:
			consensus.append(first[0])
		elif (first[1]+second[1])>log(20,2)*0.75:
			consensus.append(first[0]+'/'+second[0])
		else:
			consensus.append('*')

		aaDict={}
		height={}
	return consensus
