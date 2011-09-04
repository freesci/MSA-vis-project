from django.db import models
from django.forms import ModelForm, Textarea
from django import forms

# Create your models here.

def check_correctness(inputsequences,fields_name):
  try:
    from Bio import AlignIO
    from Bio import SeqIO
    import StringIO
    x = StringIO.StringIO(inputsequences)
    s = SeqIO.parse(x, "fasta")
  except:
    msg = "Input is not in FASTA format!"
    self._errors[fields_name] = self.error_class([msg])
  try:
    seqDict = SeqIO.to_dict(s)
  except ValueError:
    x.close()
    msg = "Input contains repeated names of sequences!"
    self._errors[fields_name] = self.error_class([msg])
  print "seqDict",seqDict
  if len(seqDict) < 2:
    x.close()
    self._errors[fields_name] = self.error_class(["Input contains less than 2 sequences!"])
  n=0
  for i in xrange(len(seqDict.values()[0].seq)):
    for key in seqDict.keys():
      if seqDict[key][i]!="-":
	n+=1
  if n==0:
    x.close()
    msg = "Wrong MSA. There is no non-gap letter on position %d!" % i
    self._errors[fields_name] = self.error_class([msg])

gchoices= (
	('Slow', 'Run PSI-BLAST'),
        ('Fast', 'without PSI-BLAST'),
        )
                
class Page(models.Model):
  sequences = models.TextField(max_length = 400,blank= True,help_text="Enter a sequence alignment in FASTA format (max sequences length <= 400)") # help_text nie jest najladniejszym rozwiazaniem w widoku formularzu...
  upload_file = models.FileField(upload_to = "uploaded_files", blank= True, help_text="Or upload a file") # uploaded_files jest katalogiem tworzonym wzgledem MEDIA_ROOT
  email = models.EmailField(max_length = 20,blank= True,help_text="Send visualization to (optional):") #blank=True, aby to pole moglo pozostac puste przy wypelnianiu pol
  unixtime = models.IntegerField(max_length = 40)
  choice = models.CharField(max_length=4, choice=gchoices, default = "Slow")
  linewidth = models.IntegerField(max_length = 300,default=30,blank=True,null = True,help_text="Number of aminoacids in one row in graph (default=30)")

class PageForm(forms.ModelForm):
  choice = forms.CharField(max_length=4, widget=forms.Select(choice=gchoices),help_text="Predict secondary structure (PSIPRED) with PSI-BLAST or without")
  
  def clean(self):
    cleaned_data = self.cleaned_data

    sequences = cleaned_data["sequences"]
    linewidth = cleaned_data["linewidth"]
    upload_file = cleaned_data["upload_file"]

    if sequences=="" and upload_file==None:
      msg = 'You must specify an input source!'
      self._errors["sequences"] = self.error_class([msg])
      del cleaned_data["sequences"]
    if sequences!="" and upload_file!=None:
      msg = 'Select alignment or upload a local file'
      self._errors["sequences"] = self.error_class([msg])
      del cleaned_data["sequences"]

    if sequences!="":
      check_correctness(sequences,"sequences")
	
    if upload_file!=None:
      content = upload_file.read()
      check_correctness(content,"upload_file")
       
    return self.cleaned_data


  class Meta:
    model = Page
    fields = ('sequences',"upload_file", 'email',"linewidth","choice") # kolejnosc wyswietlania na stronie
   