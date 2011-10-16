from msa_vis.msa_vis_app.models import Page,PageFormfrom django.shortcuts import render_to_responsefrom django.http import HttpResponseRedirect, HttpResponsefrom django.core.urlresolvers import reverseimport msa_vis.settings as settingsfrom os import pathfrom time import localtime,mktimefrom datetime import datetimeimport subprocess as subfrom msa_vis.send_mail import send_email_tosettingsMEDIA_ROOT,settingsPAGE_ADRESS = settings.MEDIA_ROOT,settings.PAGE_ADRESSdef create_picture(choice,settingsMEDIA_ROOT,file_absolutepath,jobID,linewidth):  if choice=="Slow":    return sub.Popen([settingsMEDIA_ROOT + "uploaded_files/msavisproject.py","SLOW", file_absolutepath,jobID,linewidth,settingsMEDIA_ROOT],cwd=settingsMEDIA_ROOT + "uploaded_files/results/")  elif choice=="Fast":    return sub.Popen([settingsMEDIA_ROOT + "uploaded_files/msavisproject.py","FAST", file_absolutepath,jobID,linewidth,settingsMEDIA_ROOT],cwd=settingsMEDIA_ROOT + "uploaded_files/results/")def first_page(request):    if request.method == 'POST':      now = datetime.now()      utime = int(mktime(now.timetuple()))      page = Page(unixtime = utime)      form = PageForm(request.POST,request.FILES,instance=page)      if form.is_valid():	mail = form.cleaned_data["email"]	sequences = form.cleaned_data["sequences"]	linewidth = form.cleaned_data["linewidth"]	choice = form.cleaned_data["choice"]	m = form.save()	ID = m.id	jobID = str(ID) + "-" + str(utime)	date = datetime.fromtimestamp(mktime(localtime(int(utime))))	if m.upload_file:	  file_absolutepath = m.upload_file.path	  create_picture(choice,settingsMEDIA_ROOT,file_absolutepath,jobID,str(linewidth))	  if mail !="": send_email_to(mail,settingsPAGE_ADRESS,jobID,date)	if sequences!="":	  file_absolutepath = settingsMEDIA_ROOT + "uploaded_files/filename.fasta"	  from Bio import AlignIO,SeqIO	  import StringIO	  x=StringIO.StringIO(sequences)	  output_handle = open(file_absolutepath, "w")	  SeqIO.write(SeqIO.parse(x, "fasta"), output_handle, "fasta")	  output_handle.close()	  x.close()	  create_picture(choice,settingsMEDIA_ROOT,file_absolutepath,jobID,str(linewidth))	  if mail !="": send_email_to(mail,settingsPAGE_ADRESS,jobID,date)        return HttpResponseRedirect(reverse('msa_vis.msa_vis_app.views.second_page',args=[ID,utime]))      else:	return render_to_response("first_page.html",{"form":form})    else:      form = PageForm() # create empty form    return render_to_response("first_page.html",{"form":form})def second_page(request,ID,unixtime):  jobID = str(ID) + "-" + str(unixtime)  try:    #checking from database if this id exists and unixtime's value with this id    utime = Page.objects.get(pk=ID).unixtime    #checking if unixtime from query and unixtime from database (with appropriate id) are the same    if int(unixtime)!=int(utime):      return render_to_response("notfound_page.html",{"jobID":jobID})  except Page.DoesNotExist:    return render_to_response("notfound_page.html",{"jobID":jobID})  structtime = localtime(int(utime))  date = datetime.fromtimestamp(mktime(structtime))  if not path.exists(settingsMEDIA_ROOT + "uploaded_files/results/" + "finalMSAvis" + jobID + ".svg"):    return render_to_response("between_page.html",{"jobID":jobID, "date":date})  return render_to_response("second_page.html",{"jobID":jobID,"date":date})