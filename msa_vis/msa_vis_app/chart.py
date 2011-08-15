#!/usr/bin/python
# -*- coding: utf-8 -*-

import django
from django.http import HttpResponse

#def simple(request):
  #return HttpResponse("Hello world")
  
def simple(request):
  import matplotlib
  import matplotlib.pyplot as plt
  import numpy as np
  from pylab import * # dla dzialania arange
  from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
  import datetime
  
  
  fig = plt.figure(figsize =(7,7))
  plt.axes([0.2, 0.2, 0.7, 0.6])   # leave room below the axes for the table


  dict = {'A': -4, 'B': 4, 'C': 0.1, 'D': -2.5, 'E': 3}
  def wartosc_bezwgledna(x):
    if x < 0:	return -x
    else:		return x
  maximum = max([wartosc_bezwgledna(x) for x in dict.values()])


  plt.grid(True)    #pokaz siatke
  colLabels = ('%d' % x for x in (1,2,3,4,5))
  rowLabels = ['%d year' % x for x in (100, 50, 20, 10, 5)]


  cellText = []
  yoff = array([0.0] * len(dict.values())) # the bottom values for stacked bar chart
  #print(yoff)
  #yoff=[0.0]
  i=0
  print(dict.values())
  v=[-4, 4, 0.1,-2.5, 3] #niestety w slowniku wartosci sie automatycznie sortuja..
  for value in v:
    col = ((value+maximum)/10.0, 0.3, 0.4)    #kolor slupka zalezny od wartosci
    bar(i, value,width=1,bottom=0, color=col);
    yoff = yoff + value
    cellText.append(['%1.1f' % (x/1000.0) for x in yoff])
    i+=1
    
  print(cellText)
  plt.xticks(np.arange(len(dict))+0.4, dict.keys());
  vals = arange(-maximum-1, maximum+2, 1)
  plt.yticks(vals, ['%d' % val for val in vals]);

  #opisanie wykresu
  x = datetime.datetime.now()
  dt = x.strftime("%d/%m/%y %H:%M:%S")
  
  plt.title("MSA   "+ dt)
  plt.xlabel("Amino Acid")
  plt.ylabel("Hydrophobicity")

  canvas=FigureCanvas(fig)
  response=django.http.HttpResponse(content_type='image/png')
  canvas.print_png(response)
  return response
