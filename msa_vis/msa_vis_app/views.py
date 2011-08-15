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
