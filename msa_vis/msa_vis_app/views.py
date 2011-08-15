# Create your views here.
from msa_vis.msa_vis_app.models import Page,PageForm
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse


from django.http import HttpResponse
from django.template import Context, loader


#def first_page(request):
  #try:
    ##page = Page.objects.get(pk=page_name)	`
    #sequences = page.sequences
    #upload_file = page.upload_file
    #seqID = page.seqID
  #except Page.DoesNotExist:
    #return render_to_response("first.html",{"sequences":sequences,"seqID":seqID})
  #email = page.email # pole email nie musi byc wypelnione
  #print(page_name)
  #return render_to_response("first.html",{"page_name":page_name,"content":content})
  
  
#def first_page(request):
    #latest_poll_list = Poll.objects.all().order_by('-pub_date')[:5]
    #t = loader.get_template('polls/index.html')
    #c = Context({
        #'latest_poll_list': latest_poll_list,
    #})
    #return HttpResponse(t.render(c))
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
         return HttpResponseRedirect('msa_vis/msa_vis/edit/') # Redirect after POST
      #else: # napisac wyjatek jesli uztkownik zle wypelni jakies pole
      
    else: # czyli tu tylko ogladanie strony first_page.html
      page_form = PageForm() # tu tworze pusty formularz
      return render_to_response("first_page.html",{"page_form":page_form})
    
    
    #def contact(request):
    #if request.method == 'POST': # If the form has been submitted...
        #form = ContactForm(request.POST) # A form bound to the POST data
        #if form.is_valid(): # All validation rules pass
            ## Process the data in form.cleaned_data
            ## ...
            #return HttpResponseRedirect('/thanks/') # Redirect after POST
    #else:
        #form = ContactForm() # An unbound form

    #return render_to_response('contact.html', {
        #'form': form,
    #})
