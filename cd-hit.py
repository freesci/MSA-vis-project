#!/usr/bin/env python
import os
import subprocess
from Bio import SeqIO
import re
import sys


def cd_hit(list_of_seq):
	
	#if there are less than 15 seqs it is unnecessary to use cd-hit
	if len(list_of_seq) < 15:
		return list_of_seq
	
	
	#creating current_seq_temp (necessary for cd-hit)
	current_seq_temp = open("current_seq_temp.fasta","w+")
	try:
		SeqIO.write( list_of_seq , current_seq_temp , "fasta")
	finally:
		current_seq_temp.close()
	new_length = len(list_of_seq)
	
	#list of tuples (sequence identity threshold * 100, word_length) - see more in CD-HIT Users Guide
	th = [(100,5),(95,5),(90,5),(85,5),(80,5),(75,5),(70,5),(65,4),(60,4),(55,3),(50,3),(45,2),(40,2)]

	nr = 2
	count = 0
	while (new_length < 10 or new_length > 15) and count < 50:
		count +=1
		
		if nr%2 == 0:
			ip = "current_seq_temp.fasta"
			op = "new_seq_temp.fasta"
		else:
			ip = "new_seq_temp.fasta"
			op = "current_seq_temp.fasta"
		
		# calling cd-hit
		try:
			b = subprocess.check_output(["cd-hit","-i",ip,"-o",op,"-c",str(th[nr-1][0]/100.),"-n",str(th[nr-1][1])])
		except Exception:
			print "Please install cd-hit"
			exit(1)
		# parsing output
		m = re.search('(finished).*(clusters)',b)
		
		# new_length is number of seqs which we already have (in first step it is length of list_of_seq) 
		new_length = int(m.group(0).split()[1])
		
		
		if new_length > 15:
			nr+=1
		
		# making th bigger
		elif new_length < 10:
			th.insert(nr-1,((th[nr-2][0]+th[nr-1][0])/2, ((th[nr-2][1]+th[nr-1][1])/2)%1))
		
		count+=1
		
		# if the number of sequences is still bigger than 15 we can't do anything else (see more in CD-HIT Users Guide)
		if nr == 14:
			print "There is a very little similarity between your sequences"
			break
		
	record = list(SeqIO.parse("current_seq_temp.fasta","fasta"))
	
	# removing temporary files
	
	os.remove("current_seq_temp.fasta")
	os.remove("current_seq_temp.fasta.bak.clstr")
	os.remove("current_seq_temp.fasta.clstr")
	
	os.remove("new_seq_temp.fasta")
	os.remove("new_seq_temp.fasta.bak.clstr")
	os.remove("new_seq_temp.fasta.clstr")
	
	
	return record[:15]
