#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess as sub, subprocess
import re


handle = open(sys.argv[1], "rU")
seqDict={}
for record in SeqIO.parse(handle, "fasta") :
	seqDict[record.id]=record.seq.data
handle.close()


def PredSecondaryStructure(seqDict):
  
  """ wersja v0.0 nie korzystająca z zewnętrznej białkowej bazy danych.
      Program wykorzystuje prawdopodobieństwa wyliczone przez runpsipred_simple dla wystąpienia petli,
      helisy i kartki dla odpowiednich pozycji aminokwasow w podanych sekwencjach bialkowych,
      zwraca listę stringów: "c"=coil, "h"=helix, "s"=strand """
  
  id = seqDict.keys()
  sequences = seqDict.values()

  dic={}
  for i in xrange(len(sequences[1])):
    dic[i+1]=[]
    
  for number_seq,seq in enumerate(sequences, start=0):
    where_gaps=[]  
    for t,i in enumerate(seq, start=1):
      if i=="-":
	where_gaps.append(t)

    # dla kazdej sekwencji tworzę plik w formacie fasta i wykorzystuje go przy uruchomieniu runpsipred_simple
    sequences = SeqRecord(Seq(seq), id = id[number_seq])
    tmp_filename = "tmp"+str(number_seq)+".fasta"
    output_handle = open(tmp_filename, "w")
    SeqIO.write(sequences, output_handle, "fasta")
    output_handle.close()

    p = sub.Popen(["./runpsipred_single", tmp_filename], stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.STDOUT) # psipred sam sie pozbywa gapow
    child_output, child_error = p.communicate(input="234")

    # wczytuję jeden z wynikowych plików programu psipred_single
    match = re.match("([^.]).[^.]*", tmp_filename)
    name_tmp= match.group()+".ss"
    file_tmp=open(name_tmp, "r")

    k=1 #pozycja w sekwencji z gapami
    for l,line in enumerate(file_tmp, start=1):
      pred,coil,helix,strand = line[7],float(line[11:16]),float(line[18:23]),float(line[25:30]) #przewidywana strukt. 2-go rzed.,prawd.petli,prawd.helix,prawd.kartki
      while k in where_gaps:
	dic[k].append((0,0,0)) #gap
	k+=1
      dic[k].append((coil,helix,strand))
      k+=1
    while k in where_gaps:
      dic[k].append((0,0,0)) #gap
      k+=1
  
    retcode = subprocess.call(["rm","tmp"+str(number_seq)+".fasta"])
    retcode = subprocess.call(["rm", "tmp"+str(number_seq)+".ss"])
    retcode = subprocess.call(["rm", "tmp"+str(number_seq)+".ss2"])
    retcode = subprocess.call(["rm", "tmp"+str(number_seq)+".horiz"])
  
  pred_structure=[]
  for key,value in dic.iteritems():
    coil,helix,strand=0.0,0.0,0.0
    for krotka in value:
	coil+=krotka[0]
	helix+=krotka[1]
	strand+=krotka[2]
    maximum = max(coil,helix,strand) # zakładam, że prawdop, że mogą byc takie same jest znikome
    if maximum==coil:	pred_structure.append("c")
    if maximum==helix:	pred_structure.append("h")
    if maximum==strand:	pred_structure.append("s")

  return pred_structure
