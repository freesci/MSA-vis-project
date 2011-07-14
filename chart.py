import gtk
import matplotlib
matplotlib.use('Svg')
import matplotlib.cm as cm
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
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
	wykr.set_title("MSA Visualisation", y=0.4)
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
			tmp=plt.axes([inchW, 1-(r+1)*(margin+height+seqH)-colBarH+seqH, (1-2*inchW)*(l-col*(rows-1))/float(col), height], xlabel="Amino Acid", ylabel='Hydrophobicity')
			plt.axis([(rows-1)*col-0.5, l-0.5, -5, 5])    #min and max of the x and y axes
			plt.xticks(range(col*(rows-1),l, 5))
		else:
			tmp=plt.axes([inchW, 1-(r+1)*(margin+height+seqH)-colBarH+seqH, 1-2*inchW, height], xlabel="Amino Acid", ylabel='Hydrophobicity')
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
    	win.set_title("MSA Visualisation")
        
    	canvas = FigureCanvas(fig)
    	win.add(canvas)
    	win.show_all()
    	gtk.main()
