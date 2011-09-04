#!/usr/bin/python
# -*- coding: utf-8 -*-

import textwrap
import os,sys,re
import subprocess
import subprocess as sub
from math import log

import gtk
import numpy as np

__doc__="""
Use the msa-vis-project program to visualization MSA.
You need to have a working version of
- argparse
- biopython
- termcolor
- cd-hit
- blast2, runpsipred
- gtk, matplotlib
in order to use this.

# READ ME: usuwamy powtarzajace sie sekwencje
"""

MEDIA_PATH = "/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/"


try:
	import matplotlib
except ImportError:
	print "Matplotlib packages not found. You may download it from: http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.0/"
	exit(1)

matplotlib.use('Svg')
import matplotlib.cm as cm
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas

try:
	import argparse
except ImportError:
	print "Argparse packages not found. You may download it from: http://argparse.googlecode.com/files/argparse-1.2.1.tar.gz"
	exit(1)

try:
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from Bio import SeqIO
except ImportError:
	print "Biopython packages not found. You may download it from: http://www.biosino.org/mirror/www.biopython.org/Download/default.htm"
	exit(1)

try:
	from termcolor import colored, cprint
except ImportError:
	print "Termcolor packages not found. You may download it from: http://pypi.python.org/pypi/termcolor"
	exit(1)

#Creating help
def inpout():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Visualization of MSA',prog="MSA",
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

#Checks if user gives file and if not print help 
def checkinp(parser,inp):
    if len(inp)==0:
        print colored("Error: Not specified data file\n----------------------------------",'red')
        parser.print_help()
        exit()
    return inp[0]

inp=checkinp(parser, inp)

#Creating dictionary of sequences
def seqDict(ifile):
    #Checks if in dictionary is sequence of the same name like the sequence which we wanted to add and whether they are the same
    def has_key(d,key):
        if d.has_key(key.id):
            if d[key.id]!=key.seq.data:
                key.id=key.id+"(1)"
    	return key

    #Checks if in dictionary is the same sequence like the sequence which we wanted to add
    def has_value(d,value):
       present=False
       for val in d.items():
           if val[1]==value.seq.data:
               present=True
       return present
    
    #Checks if sequence which we wanted to add is the same length as the other
    def equal_length(seq,num,lens):
        if len(seq.seq.data)!=lens:
            print colored("Error: Sequences number "+str(num)+" aren't of the same length like the other\n-----------------------------------------",'red')
            parser.print_help()
            exit()

    seqDict={}
    iteration=1
    
    #Checks if your file exists and open it
    try:
	handle = open(ifile, "rU")
    except IOError:
	print colored("Error: File '"+ifile+"' doesn't exist\n----------------------------",'red')
        exit()
    
    #Adding sequences to the dictionary
    try:
        for record in SeqIO.parse(handle, "fasta") :
            if len(seqDict)==0: 
                seqDict[record.id]=record
	        lenseq=len(record.seq.data)
            else:
		if has_value(seqDict,record)==False:
	            record=has_key(seqDict,record)
		    equal_length(record,iteration,lenseq)
	            seqDict[record.id]=record
		    lenseq=len(record.seq.data)
	    iteration+=1
    except IndexError:
        print colored("Error: Sequence without ID found\n---------------------------------------------",'red')
        parser.print_help()
        exit()

    handle.close()

    #Checks if file was empty
    if len(seqDict)==0:
        print colored("Error: The specified file is empty\n----------------------------",'red')
        parser.print_help()
        exit()     

    return seqDict

dictionary=seqDict(inp)

#Checks if number of aminoacids which user has entered is greater than the length of the sequence
if nam>len(dictionary.values()[0]):
	nam=len(dictionary.values()[0])



"""
arguments:
	list_of_seq - values of sequence dictionary
"""	
def cd_hit(dic_of_seq):

	#if there are less than 15 seqs it is unnecessary to use cd-hit
	if len(dic_of_seq) < 15:
		return dic_of_seq

	list_of_seq = dic_of_seq.values()

	for seq_record in list_of_seq:
		seq_record = seq_record.seq.ungap('-')

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
			print "Please download and install cd-hit v4.5.4 from http://code.google.com/p/cdhit/downloads/list"
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
	dic_of_rec = {}
	for x in record[:15]:
		dic_of_rec[x.name] = dic_of_seq[x.name]

	# removing temporary files

	os.remove("current_seq_temp.fasta")
	os.remove("current_seq_temp.fasta.bak.clstr")
	os.remove("current_seq_temp.fasta.clstr")

	os.remove("new_seq_temp.fasta")
	os.remove("new_seq_temp.fasta.bak.clstr")
	os.remove("new_seq_temp.fasta.clstr")


	return dic_of_rec




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
	e=19/(2*log(2)*len(seqDict[id[0]]))	#small-sample correction

	for i in xrange(len(seqDict[id[0]])):
		for key in id:
			if seqDict[key][i]!="-":
				if aaDict.has_key(seqDict[key][i]):
					aaDict[seqDict[key][i]]+=1.0
				else:
					aaDict[seqDict[key][i]]=1.0
				n+=1
		if n==0:
			print colored("Error: Wrong MSA. There is no non-gap letter on position "+str(i)+"\n-----------------------------------------",'red')
			parser.print_help()
        		exit(1)

		for k in sorted(aaDict):
			aaDict[k]/=n
			h+=aaDict[k]*log(aaDict[k], 2)

		r=log(20,2)+h-e
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

		if first[1]>(log(20,2)-e)*0.5:
			consensus.append(first[0])
		elif (first[1]+second[1])>(log(20,2)-e)*0.75:
			consensus.append(first[0]+'/'+second[0])
		else:
			consensus.append('*')

		aaDict={}
		height={}
	return consensus




""" 
    Use runpsipred (with psi-blast) or runpsipred_single (without psi-blast) program to predicting secondary structure.
    Sequence database available on http://www.ebi.ac.uk/uniprot/database/download.html.
    See more in README.
    arguments:
	seqDict - sequence dictionary
"""
def pred_secondary_structure(seqDict):
  
  id = seqDict.keys()
  sequences = seqDict.values()

  dic={}
  for i in xrange(len(sequences[0])):
    dic[i+1]=[]
   
  # runpsipred program ignores gaps in sequence, that's why it's necessary to remember theirs position in sequence
  for number_seq,seq in enumerate(sequences, start=0):
    where_no_gaps=[]  
    for t,i in enumerate(seq, start=1):
      if i!="-":	where_no_gaps.append(t)
      else:		dic[t].append((0,0,0)) # gap

    # creating current_seq_temp (necessary for runpsipred)
    output_handle = open(MEDIA_PATH + "uploaded_files/current_seq_temp.fasta", "w")
    try:
      SeqIO.write(seq, output_handle, "fasta")
    finally:
      output_handle.close()

    # calling runpsipred
    p = sub.Popen([MEDIA_PATH + "uploaded_files/runpsipred_single", MEDIA_PATH + "uploaded_files/current_seq_temp.fasta"], cwd=MEDIA_PATH + "uploaded_files/", stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.STDOUT)
    child_output, child_error = p.communicate(input="234")

    # loading one of the resulting files runpsipred
    try:
      file_tmp=open(MEDIA_PATH + "uploaded_files/current_seq_temp.ss", "r")
    except Exception:
      print "Problem with opening resulting file runpsipred."
      print "You can download and install runpsipred from http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/"
      print "If it's already installed please check variable child_error and child_output in code."
      exit(1)

    for l,line in enumerate(file_tmp, start=1):
      pred,coil,helix,strand = line[7],float(line[11:16]),float(line[18:23]),float(line[25:30])
      dic[where_no_gaps[l-1]].append((coil,helix,strand))
      
  # removing temporary files 
  l = [MEDIA_PATH + "uploaded_files/current_seq_temp.fasta",
  MEDIA_PATH + "uploaded_files/current_seq_temp.ss",
  MEDIA_PATH + "uploaded_files/current_seq_temp.ss2",
  MEDIA_PATH + "uploaded_files/current_seq_temp.horiz"]
  for i in l:  os.remove(i) 
  
  
  pred_structure=[]
  for value in dic.values():
    coil,helix,strand=0.0,0.0,0.0
    for tuple in value:
	coil+=tuple[0]
	helix+=tuple[1]
	strand+=tuple[2]
    maximum = max(coil,helix,strand)
    if maximum==coil:	pred_structure.append("c")
    if maximum==helix:	pred_structure.append("h")
    if maximum==strand:	pred_structure.append("s")

  return pred_structure





""" 
  reading file with Kyte Doolitle scale 
  
"""
def readKD():

	KD={}
	f=open(MEDIA_PATH + 'uploaded_files/KD.scale',"r")
	lines=f.readlines()
	for i in xrange(1, len(lines)):
		tmp=lines[i].split()
		KD[tmp[0]]={'v':float(tmp[3]), 'h':float(tmp[4])}

	return KD


"""
arguments:
	seqDict - sequence dictionary
	KDscale - hydrophobicity and side-chain volume for every amino acid; dictionary returned by readKD function
"""
def stat(seqDict, KDscale):
  
	id=seqDict.keys()
	hydro=[]
	chain=[]
	tmpH=0
	tmpV=0
	n=0
	for i in xrange(len(seqDict[id[0]])):
		for j in id:
			if seqDict[j][i]!='-':
				tmpH+=KDscale[seqDict[j][i]]['h']
				tmpV+=KDscale[seqDict[j][i]]['v']
				n+=1.0
		hydro.append(tmpH/n)
		chain.append(tmpV/n)
		tmpH=0
		tmpV=0
		n=0

	return hydro, chain



"""
arguments:
	consensus - consensus sequence (list)
	hydro - average hydrophobicity (list) 
	chain - average side chain volume (list)
	stru - secondary structure (list)
	cd_hit - cd-hit output (sequence dctionary)
	col - number of columns (default=30)
"""
def chart(consensus, hydro, chain, stru, cd_hit, filename, col):

	minC=60.0	#min side-chain volume
	maxC=230.0	#max side-chain volume
	l = len(consensus)
	k = len(cd_hit)
	rows = int(np.ceil(float(l)/col))	#number of rows

	fig = plt.figure(figsize =(col/3,rows*(9+k/2.0)+4))

	inchH=1.0/(rows*(9+k/2.0)+4)
	colBarH = 3.0*inchH	#colorbar height
	margin = 2.0*inchH
	height=6.0*inchH		#barchart height
	seqH = (k/2.0+1.0)*inchH	#sequence table height

	inchW = 1.0/(col/3)

	wykr=plt.axes([2*inchW, 1-colBarH, 1-4*inchW, colBarH])
	wykr.set_title("MSA Visualization", y=0.4)
	plt.axis('off')

	#side chain volume colorbar
	m=cm.ScalarMappable(cmap=cm.autumn)
	m.set_array(np.array([minC, maxC]))
	cbr=plt.colorbar(m, orientation='horizontal', fraction=0.4)
	cbr.set_label('Side Chain Volume')


	width = 1
	widths=[1.0/col]*col
	for r in xrange(rows):
		#barchart
		if  r==rows-1:
			tmp=plt.axes([2*inchW, 1-(r+1)*(margin+height+seqH)-colBarH+seqH, (1-4*inchW)*(l-col*(rows-1))/float(col), height], xlabel="Amino Acid", ylabel='Hydrophobicity')
			plt.axis([(rows-1)*col-0.5, l-0.5, -5, 5])    #min and max of the x and y axes
			plt.xticks(range(col*(rows-1),l, 5))
		else:
			tmp=plt.axes([2*inchW, 1-(r+1)*(margin+height+seqH)-colBarH+seqH, 1-4*inchW, height], xlabel="Amino Acid", ylabel='Hydrophobicity')
		    	plt.axis([r*col-0.5, (r+1)*col-0.5, -5, 5])    #min and max of the x and y axes
			plt.xticks(range(col*r,(r+1)*col, 5))

		for i in xrange(col):
			if r==rows-1 and i==l-col*(rows-1): break	#break if last chart is shorter
			c = (1, (chain[col*r+i]-minC)/(maxC-minC), 0)    #bar color
			tmp.bar(col*r+i, hydro[col*r+i], width, color=c, align='center', linewidth=1)

		if r==0:
			#consensus table
			tabCons = plt.table(cellText=[consensus[col*r:col*(r+1)]], cellLoc='center', rowLabels = ["consensus"], colWidths=widths, bbox=[0, 1.1, 1, 0.07])

			#structure table
			tabStru = plt.table(cellText=[stru[col*r:col*(r+1)]], cellLoc='center', rowLabels = ["structure"], colWidths=widths, bbox=[0, 1.02, 1, 0.07])

			#sequence table
			text=[]
			labels=[]
			for key in sorted(cd_hit.keys()):
				text.append(cd_hit[key][col*r:col*(r+1)])
				labels.append(key+": ")
			tabCdHit =  plt.table(cellText=text, cellLoc='center', colWidths=widths, rowLabels = labels, bbox=[0, -(k/2.0+1.0)/6.0, 1, k/12.0])

		elif r==rows-1:
			widths=[1.0/col]*(l-col*(rows-1))

			#consensus table
			tabCons = plt.table(cellText=[consensus[col*(rows-1):]], cellLoc='center', colWidths=widths, bbox=[0, 1.1, 1, 0.07])

			#structure table
			tabStru = plt.table(cellText=[stru[col*(rows-1):]], cellLoc='center', colWidths=widths, bbox=[0, 1.02, 1, 0.07])

			#sequence table			
			text=[]
			for key in sorted(cd_hit.keys()):
				text.append(cd_hit[key][col*r:col*(r+1)])
			tabCdHit =  plt.table(cellText=text, cellLoc='center', colWidths=widths, bbox=[0, -(k/2.0+1.0)/6.0, 1, k/12.0])
		else:
			#consensus table
			tabCons = plt.table(cellText=[consensus[col*r:col*(r+1)]], cellLoc='center', colWidths=widths, bbox=[0, 1.1, 1, 0.07])

			#structure table
			tabStru = plt.table(cellText=[stru[col*r:col*(r+1)]], cellLoc='center', colWidths=widths, bbox=[0, 1.02, 1, 0.07])

			#sequence table			
			text=[]
			for key in sorted(cd_hit.keys()):
				text.append(cd_hit[key][col*r:col*(r+1)])
			tabCdHit =  plt.table(cellText=text, cellLoc='center', colWidths=widths, bbox=[0, -(k/2.0+1.0)/6.0, 1, k/12.0])

		tabCons.auto_set_font_size(False)
		tabCons.set_fontsize(8)
		tabStru.auto_set_font_size(False)
		tabStru.set_fontsize(8)
		tabCdHit.auto_set_font_size(False)
		tabCdHit.set_fontsize(12)
		for v in tabCdHit.get_celld().values():
			v.set_edgecolor('w')

	plt.savefig(filename)
	return fig

def window(fig):
	win = gtk.Window()
	win.connect("destroy", gtk.main_quit)
	win.set_default_size(1200, 1200)
    	win.set_title("MSA Visualization")
        
    	canvas = FigureCanvas(fig)
    	win.add(canvas)
    	win.show_all()
    	gtk.main()




if __name__ == "__main__":
  
 consensus = consensus(dictionary)
 stru = pred_secondary_structure(dictionary)
 hydro, chain = stat(dictionary,readKD())
 #cd_hit = cd_hit(dictionary)
 cd_hit = dictionary

 fig = chart(consensus, hydro, chain, stru, cd_hit, out, nam)