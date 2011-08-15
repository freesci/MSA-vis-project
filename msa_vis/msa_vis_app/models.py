from django.db import models

#from django.core.files.storage import FileSystemStorage
#fs = FileSystemStorage(location='/home/kasia/Pulpit/kasia/KASIA/IBB_praktyki/MSA-vis-project/MSA-vis-project') # w koncu nei uzylam tego w argumecie funkcji models.FileField, nie wiem czy slusznie

# Create your models here.


class Page(models.Model):
  sequences = models.TextField(max_length = "150",help_text="(Max dlugosc seq <= 150)") # help_text nie jest najladniejszym rozwiazaniem w widoku formularzu...
  upload_file = models.FileField(upload_to = "uploaded_files", help_text="(Wybierz plik w formacie fasta)") #storage=fs # uploaded_files jest katalogiem tworzanym wzgledem MEDIA_ROOT
  email = models.EmailField(blank=True, max_length = "20",help_text="(narazie jest obowiazkowe)") #blank=True, aby to pole moglo pozostac puste przy wypelnianiu pol
  #seqID = models.CharField(max_length = "40") # podobe pole do tego ze strony aln2plot - tam jesli uzytkownik nic w to pole nie wipisze, program sam wygeneruje losowe id
  
  def __unicode__(self): # ?!
    return self.sequences
#asia@kasia-laptop:~/Pulpit/kasia/KASIA/IBB_praktyki/django/msa_vis$ ./manage.py syncdb 
#Creating table msa_vis_app_page


from django.forms import ModelForm

class PageForm(ModelForm):
  class Meta:
    model = Page # w klasie meta mus byc pole o nazwie model, do jakiego modelu ma byc stworzony formularz
    fields = ('sequences', 'upload_file','email') # w tej kolejnosci pokaza sie te rzeczy w formularzu na stronie
    
    
    