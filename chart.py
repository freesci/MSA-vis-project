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
	col - number of columns (default=30)
"""

def chart(consensus, hydro, chain, stru, col=50):
	minC=60.0	#min side-chain volume
	maxC=230.0	#max side-chain volume
	l = len(consensus)
	if col==0: col=l
	rows = int(np.ceil(float(l)/col))+1	#number of rows
	margin = 0.1
	height=(1-margin*(rows+1))/rows		#barchart height

	fig = plt.figure(figsize =(col/3,rows*8))

	wykr=plt.axes([margin, 1-(margin+height), 0.8, height])
	wykr.set_title("MSA Visualisation", y=1)
	plt.axis('off')

	#side chain volume colorbar
	m=cm.ScalarMappable(cmap=cm.autumn)
	m.set_array(np.array([minC, maxC]))
	cbr=plt.colorbar(m, orientation='horizontal', fraction=0.7)
	cbr.set_label('Side Chain Volume')


	width = 1
	widths=[1.0/col]*col
	for r in xrange(rows-2):

		#barchart
		tmp=plt.axes([margin, 1-(r+2)*(margin+height), 1-2*margin, height], xlabel="Amino Acid", ylabel='Hydrophobicity')
		for i in xrange(col):
			c = (1, (chain[col*r+i]-minC)/(maxC-minC), 0)    #bar color
			tmp.bar(col*r+i, hydro[col*r+i], width, color=c, align='center', linewidth=1)
	    	plt.axis([r*col-0.5, (r+1)*col-0.5, -5, 5])    #min and max of the x and y axes
		if r==0:
			#consensus table
			tabCons = plt.table(cellText=[consensus[col*r:col*(r+1)]], cellLoc='center', rowLabels = ["consensus"], colWidths=widths, bbox=[0, 1.1, 1, 0.07])
			tabCons.auto_set_font_size(False)
			tabCons.set_fontsize(8)

			#structure table
			tabStru = plt.table(cellText=[stru[col*r:col*(r+1)]], cellLoc='center', rowLabels = ["structure"], colWidths=widths, bbox=[0, 1.02, 1, 0.07])
			tabStru.auto_set_font_size(False)
			tabStru.set_fontsize(8)
		else:
			#consensus table
			tabCons = plt.table(cellText=[consensus[col*r:col*(r+1)]], cellLoc='center', colWidths=widths, bbox=[0, 1.1, 1, 0.07])
			tabCons.set_fontsize(8)

			#structure table
			tabStru = plt.table(cellText=[stru[col*r:col*(r+1)]], cellLoc='center', colWidths=widths, bbox=[0, 1.02, 1, 0.07])
			tabStru.set_fontsize(8)

	#last (shorter) barchart
	tmp=plt.axes([margin, margin, (1-2*margin)*(l-col*(rows-2))/float(col), height], xlabel="Amino Acid", ylabel='Hydrophobicity')
	for i in xrange(l-col*(rows-2)):
		c = (1, (chain[col*(rows-2)+i]-minC)/(maxC-minC), 0)    #bar color
		tmp.bar(col*(rows-2)+i, hydro[col*(rows-2)+i], width, color=c, align='center', linewidth=1)
	plt.axis([(rows-2)*col-0.5, l-0.5, -5, 5])    #min and max of the x and y axes
	plt.xticks(range(col*(rows-2),l, 5))
    	
	widths=[1.0/col]*(l-col*(rows-2))
	tabCons = plt.table(cellText=[consensus[col*(rows-2):]], cellLoc='center', colWidths=widths, bbox=[0, 1.1, 1, 0.07])
	tabCons.auto_set_font_size(False)
	tabCons.set_fontsize(8)

	tabStru = plt.table(cellText=[stru[col*(rows-2):]], cellLoc='center', colWidths=widths, bbox=[0, 1.02, 1, 0.07])
	tabStru.auto_set_font_size(False)
	tabStru.set_fontsize(8)
	
	plt.savefig('test')
	return fig

def window(fig):
	win = gtk.Window()
	win.connect("destroy", gtk.main_quit)
	win.set_default_size(1200, 1200)
    	win.set_title("MSA Visualisation")
        
    	canvas = FigureCanvas(fig)
    	win.add(canvas)
    	win.show_all()
    	gtk.main()
