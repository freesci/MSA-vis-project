from django.db import models
from django.forms import ModelForm

# Create your models here.

class Page(models.Model):
  sequences = models.TextField(max_length = "150",blank= True,help_text="(Max dlugosc seq <= 150)") # help_text nie jest najladniejszym rozwiazaniem w widoku formularzu...
  upload_file = models.FileField(upload_to = "uploaded_files", blank= True, help_text="(Wybierz plik w formacie fasta)") # uploaded_files jest katalogiem tworzonym wzgledem MEDIA_ROOT
  email = models.EmailField(max_length = "20",blank= True,help_text="(pole nieobowiazkowe)") #blank=True, aby to pole moglo pozostac puste przy wypelnianiu pol
  seqID = models.CharField(max_length = "40",blank= True,help_text="(pole nieobowiazkowe)") # podobe pole do tego ze strony aln2plot - tam jesli uzytkownik nic w to pole nie wipisze, program sam wygeneruje losowe id

class PageForm(ModelForm):
  #def clean_sequences(self):        # walidacja na serwerze? ....
    #data = self.cleaned_data['sequences']
    #if "cos" not in data:
      #raise forms.ValidationError("You have forgotten about cos!")
  class Meta:
    model = Page # w klasie meta mus byc pole o nazwie model, do jakiego modelu ma byc stworzony formularz
    fields = ('sequences', 'email','seqID',"upload_file") # kolejnosc wyswietlania na stronie
    
    