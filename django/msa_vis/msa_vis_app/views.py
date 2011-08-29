# Create your views here.
from msa_vis.msa_vis_app.models import Page,PageForm
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.template import Context, loader
from models import PageForm
import msa_vis.settings as settings


def check_correctnessfile(absolutefilepath):  # spr czy wprowadzone msa jest prawidlowe
	  try:
	      from Bio.Seq import Seq
	      from Bio.SeqRecord import SeqRecord
	      from Bio import SeqIO
	      #http://biopython.org/DIST/docs/tutorial/Tutorial.html#SeqIO:to_dict
	      try:
	      	seqDict = SeqIO.to_dict(SeqIO.parse(absolutefilepath, "fasta")) #spr czy zaladowany plik jest w formacie .fasta
	      except ValueError:
	      	#error1 = "Input contains repeated names of sequences!"
	      	return False
	      n=0
	      if len(seqDict) < 2: 
		#error2 = "Input contains less than 2 sequences!"
		return False
	      for i in xrange(len(seqDict.values()[0].seq)):
		for key in seqDict.keys():
		  if seqDict[key][i]!="-":
		    n+=1
	      if n==0:
	      	#error3 =  "Wrong MSA. There is no non-gap letter on position" + i
	      	return False 
	      return True
	  except:
	  	#error4 = "Wrong MSA"
	  	return False

def first_page(request):
   if request.method == 'POST':
     form = PageForm(request.POST,request.FILES)
     if form.is_valid():
      
	mail = form.cleaned_data["email"]
	sequences = form.cleaned_data["sequences"]
	seqID = form.cleaned_data["seqID"]
	linewidth = form.cleaned_data["linewidth"]

	m=form.save()
	if seqID=="":
	  import random
	  seqID =  random.randint(1, 100000)
	if linewidth is not None:  linewidth= str(linewidth)
	# mozliwe uzycie albo upload fieldu albo sequences jednoczesnie
	if m.upload_file and sequences:		return HttpResponseRedirect('error/')
	if m.upload_file: # tu problemem moze byc to, ze jesli uzytkownik zaladuje mnostwo razy plik o tej samej nazwie (wtedy w bazie dodawane sa kolejne "_" do nazwy pliku - ich dodawanie jest pewnie ograniczone) -  to moze sie zapchac baza..	  if not check_correctnessfile(m.upload_file.path): return HttpResponseRedirect('error1/')
	  file_absolutepath = m.upload_file.path
	  if not check_correctnessfile(file_absolutepath): return HttpResponseRedirect('error/')
	  import subprocess as sub #rozwiazanie na teraz; latwo zapchac komputer jesli uztkownik odpali wiele procesow; system kolejkowania - zajrzec w google: celery
	  if linewidth is not None:
	    p = sub.Popen([settings.MEDIA_ROOT+"uploaded_files/msavisproject.py", file_absolutepath, "-o","MSAvis"+str(seqID)+".svg", "-a",linewidth],cwd=settings.MEDIA_ROOT+"uploaded_files/")
	  else:
	    p = sub.Popen([settings.MEDIA_ROOT+"uploaded_files/msavisproject.py", file_absolutepath, "-o","MSAvis"+str(seqID)+".svg"],cwd=settings.MEDIA_ROOT+"uploaded_files/")
	  child_output, child_error = p.communicate(input="234")
	else: # jesli nie zostal zaladowany plik w UploadField, to oczekuje, ze zostala wprowadzana sekwencja do field'u sequences
	
	  if sequences!="":
	    f = settings.MEDIA_ROOT+"uploaded_files/filename.fasta"
	    try:
	      from Bio import AlignIO
	      from Bio import SeqIO
	      import StringIO
	      x=StringIO.StringIO(sequences)
	      output_handle = open(f, "w")
	      s = SeqIO.parse(x, "fasta")
	      count = SeqIO.write(s, output_handle, "fasta")
	      output_handle.close()
	      x.close()
	    except:
	      return HttpResponseRedirect('error/')
	    if not check_correctnessfile(f): return HttpResponseRedirect('error/')
	    import subprocess as sub #rozwiazanie na teraz; latwo zapchac komputer jesli uztkownik odpali wiele procesow; system kolejkowania - zajrzec w google: celery
	    if linewidth is not None:
	      p = sub.Popen([settings.MEDIA_ROOT+"uploaded_files/msavisproject.py", f, "-o","MSAvis"+str(seqID)+".svg","-a",linewidth],cwd=settings.MEDIA_ROOT+"uploaded_files/")
	    else:
	      p = sub.Popen([settings.MEDIA_ROOT+"uploaded_files/msavisproject.py", f, "-o","MSAvis"+str(seqID)+".svg"],cwd=settings.MEDIA_ROOT+"uploaded_files/")
	    child_output, child_error = p.communicate(input="234")
	  else:
	    return HttpResponseRedirect('error/')
	 
	#gdyby stara baza danych uniemozliwiala pojscie dalej - np. Type error i jakies krzaki
	#/manage.py reset msa_vis_app
	# usuniecie wszsytkich tabel dotyczacych aplikacji z bazy danych i utworzenie ich na nowo
	# ./manage syncdb
		

	if mail!="":
	  from django.core.mail import EmailMessage
	  text_content = 'This is an important message.'
	  html_content = '<p>This is an <strong>important</strong> message.</p>'
	  msg = EmailMessage('subject',html_content,to = [mail])
	  msg.content_subtype = "html"
	  msg.attach_file(settings.MEDIA_ROOT+"uploaded_files/MSAvis"+str(seqID)+".svg")
	  msg.send()

        return render_to_response("second_page.html",{"seqID":seqID})
     else: # wyjatek jesli uztkownik zle wypelni jakies pole
	return HttpResponseRedirect('error/')
   else: # tu tylko ogladanie strony first_page.html
      page_form = PageForm() # tu tworze pusty formularz
      return render_to_response("first_page.html",{"page_form":page_form})
 
def second_page(request,seqID):
  return render_to_response("second_page.html",{"seqID":seqID})
  
def third_page(request):
  return render_to_response("third_page.html",{})