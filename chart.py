import gtk
import random
import matplotlib
matplotlib.use('Svg')
import matplotlib.cm as cm
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.transforms as trans
import numpy as np
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas



def readKD():
	KD={}
	f=open('KD.scale')
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
"""

def chart(consensus, hydro, chain, stru):
	col=30	#aas in row
	minC=60.0	#min side-chain volume
	maxC=230.0	#max side-chain volume
	l = len(consensus)
	rows = l/col

	fig = plt.figure(figsize =(col/3,(rows+1)*8))
	plt.subplots_adjust(hspace=0.7)


	wykr=fig.add_subplot(rows+2, 1, 1, axisbg=None)
	wykr.set_title("MSA Visualisation", y=1)

	#side chain volume colorbar
	m=cm.ScalarMappable(cmap=cm.autumn)
	m.set_array(np.array([minC, maxC]))
	cbr=plt.colorbar(m, orientation='horizontal', fraction=0.7)
	cbr.set_label('Side Chain Volume')
	plt.axis('off')	

	width = 1
	widths=[1.0/col]*col
	for r in xrange(rows):

		#barchart
		tmp=fig.add_subplot(rows+2, 1, r+2, xlabel="Amino Acid", ylabel='Hydrophobicity')
		for i in xrange(col):
			#print "r, i", r, i
			c = (1, (chain[col*r+i]-minC)/(maxC-minC), 0)    #bar color
			tmp.bar(col*r+i, hydro[col*r+i], width, color=c, align='center', linewidth=1)
	    	plt.axis([r*col-0.5, (r+1)*col-0.5, -5, 5])    #min and max of the x and y axes
		
		#consensus table
		rnames=["consensus"]
		celltext=[consensus[col*r:col*(r+1)]]
		tabCons = plt.table(cellText=celltext, cellLoc='center', rowLabels = rnames, colWidths=widths, bbox=[0, 1.1, 1, 0.07])
		tabCons.set_fontsize(8)

		#structure table
		rnames=["structure"]
		celltext=[stru[col*r:col*(r+1)]]
		tabStru = plt.table(cellText=celltext, cellLoc='center', rowLabels = rnames, colWidths=widths, bbox=[0, 1.02, 1, 0.07])
		tabStru.set_fontsize(8)

	#TODO: w oknie wyswietla sie dobrze, ale w zapisanym obrazku tabele sa rozciagniate i przesuniete + cos sie psuje z rozmiarem czcionki
	if float(col)/(l-col*rows)>2:
		cols=np.ceil(float(col)/(l-col*rows))
	else:
		cols=1
	
	tmp=fig.add_subplot(rows+2, cols, (rows+1)*cols+1, xlabel="Amino Acid", ylabel='Hydrophobicity')
	for i in xrange(l-col*rows):
		#print "r, i", rows, i
		c = (1, (chain[col*rows+i]-minC)/(maxC-minC), 0)    #bar color
		tmp.bar(col*rows+i, hydro[col*rows+i], width, color=c, align='center', linewidth=1)
	plt.axis([rows*col-0.5, l-0.5, -5, 5])    #min and max of the x and y axes
	plt.xticks(range(col*rows,l, 5))
    	
	rnames=["consensus"]
	celltext=[consensus[col*rows:]]
	widths=[1.0/col]*(l-col*rows)
	tabCons = plt.table(cellText=celltext, cellLoc='center', rowLabels = rnames, colWidths=widths, bbox=[0, 1.1, 1, 0.07])
	#tabCons.set_fontsize(8)

	rnames=["structure"]
	celltext=[stru[col*rows:]]
	tabStru = plt.table(cellText=celltext, cellLoc='center', rowLabels = rnames, colWidths=widths, bbox=[0, 1.02, 1, 0.07])
	#tabStru.set_fontsize(8)
	plt.savefig('test')
	return fig

def okno(fig):
	win = gtk.Window()
	win.connect("destroy", gtk.main_quit)
	win.set_default_size(800, 1200)
    	win.set_title("MSA Visualisation")
        
    	canvas = FigureCanvas(fig)
    	win.add(canvas)
    	win.show_all()
    	gtk.main()
