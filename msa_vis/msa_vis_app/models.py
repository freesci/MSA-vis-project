from django.db import models
from django.forms import ModelForm

# Create your models here.


class Page(models.Model):
  sequences = models.TextField(max_length = "150",help_text="(Max dlugosc seq <= 150)") # help_text nie jest najladniejszym rozwiazaniem w widoku formularzu...
  #upload_file = models.FileField(upload_to = "uploaded_files", help_text="(Wybierz plik w formacie fasta)") #storage=fs # uploaded_files jest katalogiem tworzanym wzgledem MEDIA_ROOT
  email = models.EmailField(blank=True, max_length = "20",help_text="(chwilowo obowiazkowe)") #blank=True, aby to pole moglo pozostac puste przy wypelnianiu pol
  seqID = models.CharField(max_length = "40",help_text="(chwilowo obowiazkowe)") # podobe pole do tego ze strony aln2plot - tam jesli uzytkownik nic w to pole nie wipisze, program sam wygeneruje losowe id

class PageForm(ModelForm):
  class Meta:
    model = Page # w klasie meta mus byc pole o nazwie model, do jakiego modelu ma byc stworzony formularz
    fields = ('sequences','email','seqID') # w tej kolejnosci pokaza sie te rzeczy w formularzu na stronie
    
    
    