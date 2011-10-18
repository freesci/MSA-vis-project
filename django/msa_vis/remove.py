from django.core.management import setup_environ
import settings
from time import *
import os

setup_environ(settings)

from msa_vis_app.models import Page

def remove_img():
  #old = time()		#now
  #old = time()-86400	#one day
  old = time()-604800	#one week
  #old = time()-1209600	#two weeks
  #old = time()-2592000	#one month
  for p in Page.objects.all():
    unixtime = int(p.unixtime)
    if unixtime <= old:
      if os.path.exists(settings.MEDIA_ROOT+"uploaded_files/results/finalMSAvis"+str(p.id)+"-"+str(unixtime)+".svg"):	 
	os.remove(settings.MEDIA_ROOT+"uploaded_files/results/finalMSAvis"+str(p.id)+"-"+str(unixtime)+".svg")
      p.delete()


remove_img()