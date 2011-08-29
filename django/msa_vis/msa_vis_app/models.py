from django.db import models
from django.forms import ModelForm

# Create your models here.

class Page(models.Model):
  sequences = models.TextField(max_length = 150,blank= True,help_text="(Max dlugosc seq <= 150)") # help_text nie jest najladniejszym rozwiazaniem w widoku formularzu...
  upload_file = models.FileField(upload_to = "uploaded_files", blank= True, help_text="(Wybierz plik w formacie fasta)") # uploaded_files jest katalogiem tworzonym wzgledem MEDIA_ROOT
  email = models.EmailField(max_length = 20,blank= True,help_text="(pole nieobowiazkowe)") #blank=True, aby to pole moglo pozostac puste przy wypelnianiu pol
  unixtime = models.IntegerField(blank= True) # podobe pole do tego ze strony aln2plot - tam jesli uzytkownik nic w to pole nie wipisze, program sam wygeneruje losowe id
  linewidth = models.IntegerField(max_length = 300,default=30,blank=True,null = True,help_text="(pole nieobowiazkowe); number of aminoacids in one row in graph, default=30")

class PageForm(ModelForm):
  class Meta:
    model = Page
    fields = ('sequences', 'email',"linewidth","upload_file") # kolejnosc wyswietlania na stronie
    
    