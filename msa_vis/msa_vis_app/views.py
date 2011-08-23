<<<<<<< HEAD
# Create your views here.from msa_vis.msa_vis_app.models import Page,PageFormfrom django.shortcuts import render_to_responsefrom django.http import HttpResponseRedirectfrom django.http import HttpResponsefrom django.template import Context, loaderfrom models import PageFormdef check_correctnessfile(absolutefilepath):  # spr czy wprowadzone msa jest prawidlowe	  try:	      from Bio.Seq import Seq	      from Bio.SeqRecord import SeqRecord	      from Bio import SeqIO	      #http://biopython.org/DIST/docs/tutorial/Tutorial.html#SeqIO:to_dict	      seqDict = SeqIO.to_dict(SeqIO.parse(absolutefilepath, "fasta")) #spr czy zaladowany plik jest w formacie .fasta	      n=0	      for i in xrange(len(seqDict.values()[0].seq)):		for key in seqDict.keys():		  if seqDict[key][i]!="-": # jesli na jakiejs pozycji we wszystkich sekwencjach jest "-", msa jest nieprawidlowe		    n+=1	      if n==0:	return False # Wrong MSA. There is no non-gap letter on position i 	      return True	  except:	return Falsedef first_page(request):    if request.method == 'POST': # If the form has been submitted... # wysylanie wypelnione dane z formularza, a gdy 'GET" to to w innych sytuacjach, gdy np oglada uzytkownik strone, patrz protokol HTTP      form = PageForm(request.POST,request.FILES) # A form bound to the a data      if form.is_valid(): # All validation rules pass      	mail = form.cleaned_data["email"]	sequences = form.cleaned_data["sequences"]	seqID = form.cleaned_data["seqID"]	linewidth = form.cleaned_data["linewidth"]	m=form.save()	if linewidth is not None:  linewidth= str(linewidth)	# mozliwe uzycie albo upload fieldu albo sequences jednoczesnie	print linewidth	if m.upload_file and sequences:		return HttpResponseRedirect('error1/')	if m.upload_file: # tu problemem moze byc to, ze jesli uzytkownik zaladuje mnostwo razy plik o tej samej nazwie (wtedy w bazie dodawane sa kolejne "_" do nazwy pliku - ich dodawanie jest pewnie ograniczone) -  to moze sie zapchac baza..	  if not check_correctnessfile(m.upload_file.path): return HttpResponseRedirect('error1/')	  file_absolutepath = m.upload_file.path	  import subprocess as sub #rozwiazanie na teraz; latwo zapchac komputer jesli uztkownik odpali wiele procesow; system kolejkowania - zajrzec w google: celery	  if linewidth is not None:	    p = sub.Popen(["/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/msavisproject.py", file_absolutepath, "-o","MSAvis.svg", "-a",linewidth],cwd="/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/")	  else:	    p = sub.Popen(["/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/msavisproject.py", file_absolutepath, "-o","MSAvis.svg"],cwd="/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/")	  child_output, child_error = p.communicate(input="234")	else: # jesli nie zostal zaladowany plik w UploadField, to oczekuje, ze zostala wprowadzana sekwencja do field'u sequences		  if sequences!="":	    f = "/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/filename.fasta"	    try:	      from Bio import AlignIO	      from Bio import SeqIO	      import StringIO	      x=StringIO.StringIO(sequences) # w locie tworze sekwencje .fasta, nie zapisuje na dysku	      output_handle = open(f, "w")	      s = SeqIO.parse(x, "fasta")	      count = SeqIO.write(s, output_handle, "fasta")	      output_handle.close()	      x.close()	    except:	      return HttpResponseRedirect('error1/')	    if not check_correctnessfile(f): return HttpResponseRedirect('error1/')	    import subprocess as sub #rozwiazanie na teraz; latwo zapchac komputer jesli uztkownik odpali wiele procesow; system kolejkowania - zajrzec w google: celery	    if linewidth is not None:	      p = sub.Popen(["/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/msavisproject.py", f, "-o","MSAvis.svg","-a",linewidth],cwd="/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/")	    else:	      p = sub.Popen(["/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/msavisproject.py", f, "-o","MSAvis.svg"],cwd="/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/")	    child_output, child_error = p.communicate(input="234")	  else:	    return HttpResponseRedirect('error1/')	 	#gdyby stara baza danych uniemozliwiala pojscie dalej - np. Type error i jakies krzaki	#/manage.py reset msa_vis_app	# usuniecie wszsytkich tabel dotyczacych aplikacji z bazy danych i utworzenie ich na nowo	# ./manage syncdb			if seqID=="":	  import random	  seqID = random.random()	if mail!="":	  from django.core.mail import EmailMessage	  text_content = 'This is an important message.'	  html_content = '<p>This is an <strong>important</strong> message.</p>'	  msg = EmailMessage('subject',html_content,to = [mail])	  msg.content_subtype = "html"	  msg.attach_file("/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/MSAvis.svg") #sciezka bezwgledna do pliku	  msg.send()        return render_to_response("second_page.html",{"seqID":seqID})      else: # napisac wyjatek jesli uztkownik zle wypelni jakies pole	return HttpResponseRedirect('error/')    else: # tu tylko ogladanie strony first_page.html      page_form = PageForm() # tu tworze pusty formularz      return render_to_response("first_page.html",{"page_form":page_form}) def second_page(request,seqID):  return render_to_response("second_page.html",{"seqID":seqID})  def third_page(request):  return render_to_response("third_page.html",{})  
=======
# Create your views here.
from msa_vis.msa_vis_app.models import Page,PageForm
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse


from django.http import HttpResponse
from django.template import Context, loader


from django import forms
from models import PageForm


def check_file(relativepathfile):
	  absolutefilepath = relativepathfile.path
	  try:	#spr czy zaladowany plik jest w formacie .fasta
	    from Bio import SeqIO
	    i = 0
	    lenghtseq = 0
	    for seq_record in SeqIO.parse(file_path, "fasta"):
	      # spr czy sekw w msa maja ta sama dlugosc
	      if i==0:
		lengthseq += len(seq_record)
		i+=1
	      else:
		i+=1
		if lengthseq != len(seq_record):
		  return HttpResponseRedirect('error1/')
	  except:
	    return HttpResponseRedirect('error1/')

def first_page(request):
    if request.method == 'POST': # If the form has been submitted...
      form = PageForm(request.POST,request.FILES) # A form bound to the a data
      print form.files
      if form.is_valid(): # All validation rules pass
	mail = form.cleaned_data["email"]
	sequences = form.cleaned_data["sequences"]
	seqID = form.cleaned_data["seqID"]
	m=form.save()
	# mozliwe uzycie albo upload fieldu albo sequences jednoczesnie
	if m.upload_file and sequences:		return HttpResponseRedirect('error1/')
	if m.upload_file: # tu problemem moze byc to, ze jesli uzytkownik zaladuje mnostwo razy plik o tej samej nazwie (wtedy w bazie dodawane sa kolejne "_" do nazwy pliku - ich dodawanie jest pewnie ograniczone) -  to moze sie zapchac baza..
	  check_file(m.upload_file)
	  import os
	  filename = os.path.basename(m.upload_file.path)
	  import subprocess as sub #rozwiazanie na teraz; latwo zapchac komputer jesli uztkownik odpali wiele procesow; system kolejkowania - zajrzec w google: celery
	  p = sub.Popen(["/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/msavisproject.py", filename, "-o","msavis.svg"])
	  child_output, child_error = p.communicate(input="234")
	else: # jesli nie zostal zaladowany plik w UploadField, to oczekuje, ze zostala wprowadzana sekwencja do field'u sequences
	  if sequences!="":
	    try:
	      from Bio import AlignIO
	      from Bio import SeqIO
	      import StringIO
	      f = "/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/filename.fasta"
	      # btw znalazlam w dokumentacji http://biopython.org/DIST/docs/tutorial/Tutorial.html#SeqIO:to_dict  :P
	      x=StringIO.StringIO(sequences) # w locie tworze sekwencje .fasta, nie zapisuje na dysku
	      output_handle = open(f, "w")
	      s = SeqIO.parse(x, "fasta")
	      count = SeqIO.write(s, output_handle, "fasta")
	      output_handle.close()
	      x.close()
	      import subprocess as sub #rozwiazanie na teraz; latwo zapchac komputer jesli uztkownik odpali wiele procesow; system kolejkowania - zajrzec w google: celery
	      p = sub.Popen(["/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/msavisproject.py", "filename.fasta", "-o","msavis.svg"])
	      child_output, child_error = p.communicate(input="234")
	    except:
	      return HttpResponseRedirect('error1/')	     
	  else:
	    return HttpResponseRedirect('error1/')
	 
	#gdby stara baza danych uniemozliwiala pojscie dalej - Type error i jakies krzaki
	#/manage.py reset msa_vis_app
	# usuniecie wszsytkich tabel dotyczacych aplikacji z bazy danych i utworzenie ich na nowo
	# ./manage syncdb
		
	if seqID=="":
	  import random
	  seqID = random.random()
	if mail!="":
	  from django.core.mail import EmailMessage
	  text_content = 'This is an important message.'
	  html_content = '<p>This is an <strong>important</strong> message.</p>'
	  msg = EmailMessage('subject',html_content,to = [mail])
	  msg.content_subtype = "html"
	  msg.attach_file("/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/media/uploaded_files/filename.svg") #sciezka bezwgledna do pliku
	  msg.send()

        return render_to_response("second_page.html",{"seqID":seqID})
      else: # napisac wyjatek jesli uztkownik zle wypelni jakies pole
	return HttpResponseRedirect('error1/')
	#return render_to_response("third_page.html",{})
    else: # czyli tu tylko ogladanie strony first_page.html
      page_form = PageForm() # tu tworze pusty formularz
      return render_to_response("first_page.html",{"page_form":page_form})
 
def second_page(request,seqID):
  return render_to_response("second_page.html",{"seqID":seqID})
  
def third_page(request):
  return render_to_response("third_page.html",{})
  
>>>>>>> 26a1c9aafadb46ba8aa0d161f8246994477e31f5
