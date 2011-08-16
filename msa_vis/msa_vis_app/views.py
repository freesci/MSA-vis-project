<<<<<<< HEAD
# Create your views here.from msa_vis.msa_vis_app.models import Page,PageFormfrom django.shortcuts import render_to_responsefrom django.http import HttpResponseRedirectfrom django.core.urlresolvers import reversefrom django.http import HttpResponsefrom django.template import Context, loaderfrom django import formsfrom models import PageFormdef first_page(request):    if request.method == 'POST': # If the form has been submitted... # wysylanie wypelnione dane z formularza, a gdy 'GET" to to w innych sytuacjach, gdy np oglada uzytkownik strone, patrz protokol HTTP      form = PageForm(request.POST) # A form bound to the a data      print("/n",form.as_table)      if form.is_valid(): # All validation rules pass	email = form.cleaned_data["email"]	sequences = form.cleaned_data["sequences"]	seqID = form.cleaned_data["seqID"]	print(dir(form))	form.save()	#stara baza danych uniemozliwiala pojscie dalej - Type error i jakies krzaki	#/manage.py reset msa_vis_app	# usuniecie wszsytkich tabel dotyczacych aplikacji z bazy danych i utworzenie ich na nowo	print(dir(form))	#from Bio.Seq import Seq	#from Bio.SeqRecord import SeqRecord	#from Bio import SeqIO	#output_handle = open("seqmsadjango.fasta", "w")        #SeqIO.write(unicode(sequences), output_handle, "fasta")			#print email	#print(dir(form))	#from django.core.mail import send_mail	#send_mail('Subject here', 'Here is the message.', 'nashiraV@gmail.com',[email], fail_silently=False)	#form.save() #zapisanie w bazie danych      	#print("bk")       ## Process the data in form.cleaned_data       # powalczyc z tym file:///usr/share/doc/python-django-doc/html/ref/forms/api.html#binding-uploaded-files       ## ...        #return HttpResponseRedirect('edit/') # Redirect after POST        return render_to_response("second_page.html",{"email":email,"seqID":seqID})      else: # napisac wyjatek jesli uztkownik zle wypelni jakies pole	#print(form.cleaned_data["email"],"poelseie")	return HttpResponseRedirect('error1/')	#return render_to_response("third_page.html",{})    else: # czyli tu tylko ogladanie strony first_page.html      page_form = PageForm() # tu tworze pusty formularz      return render_to_response("first_page.html",{"page_form":page_form})#raise ValueError, "The view %s.%s didn't return an HttpResponse object." % (callback.__module__, view_name) #dopisac warunek co jeslki uzytkownik poda niepoprawne danefrom msa_vis.msa_vis_app.chart import * def second_page(request,email):  return render_to_response("second_page.html",{"email":email})  def third_page(request):  return render_to_response("third_page.html",{})  
=======
# Create your views here.
from msa_vis.msa_vis_app.models import Page,PageForm
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse


from django.http import HttpResponse
from django.template import Context, loader


from django import forms
from forms import ContactForm
def first_page(request):
    if request.method == 'POST': # If the form has been submitted... # wysylanie wypelnione dane z formularza, a gdy 'GET" to to w innych sytuacjach, gdy np oglada uzytkownik strone, patrz protokol HTTP
      print("upa")
      form = ContactForm(request.POST) # A form bound to the a data
      if form.is_valid(): # All validation rules pass
       ## Process the data in form.cleaned_data
       # powalczyc z tym file:///usr/share/doc/python-django-doc/html/ref/forms/api.html#binding-uploaded-files
       ## ...
         return HttpResponseRedirect('edit/') # Redirect after POST
      #else: # napisac wyjatek jesli uztkownik zle wypelni jakies pole
      
    else: # czyli tu tylko ogladanie strony first_page.html
      page_form = PageForm() # tu tworze pusty formularz
      return render_to_response("first_page.html",{"page_form":page_form})
      
def second_page(request):
  return render_to_response("second_page.html",{})
>>>>>>> e56b88a41b472995b36edd204e01b279ca736509
