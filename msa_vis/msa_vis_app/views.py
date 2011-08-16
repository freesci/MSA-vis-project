# Create your views here.
from msa_vis.msa_vis_app.models import Page,PageForm
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse


from django.http import HttpResponse
from django.template import Context, loader


from django import forms
from models import PageForm
def first_page(request):
    if request.method == 'POST': # If the form has been submitted... # wysylanie wypelnione dane z formularza, a gdy 'GET" to to w innych sytuacjach, gdy np oglada uzytkownik strone, patrz protokol HTTP
      form = PageForm(request.POST) # A form bound to the a data
      print("/n",form.as_table)
      if form.is_valid(): # All validation rules pass
	email = form.cleaned_data["email"]
	sequences = form.cleaned_data["sequences"]
	seqID = form.cleaned_data["seqID"]
	#print(dir(form))
	form.save()
        return render_to_response("second_page.html",{"email":email,"seqID":seqID})
      else: # napisac wyjatek jesli uztkownik zle wypelni jakies pole
	return HttpResponseRedirect('error1/')
    else: # czyli tu tylko ogladanie strony first_page.html
      page_form = PageForm() # tu tworze pusty formularz
      return render_to_response("first_page.html",{"page_form":page_form})

def second_page(request,email):
  return render_to_response("second_page.html",{"email":email})
  
def third_page(request):
  return render_to_response("third_page.html",{})
 