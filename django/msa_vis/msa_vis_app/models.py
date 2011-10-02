from django.db import models
from django.forms import ModelForm
from django import forms

# Create your models here.

gchoices= (
	('Slow', 'Slower - Run PSI-BLAST'),
        ('Fast', 'Faster - without PSI-BLAST'),
        )
                
class Page(models.Model):
  sequences = models.TextField(max_length = 400,blank= True,help_text="Enter a sequence alignment in FASTA format (max sequences length <= 400)") # help_text nie jest najladniejszym rozwiazaniem w widoku formularzu...
  upload_file = models.FileField(upload_to = "uploaded_files", blank= True, help_text="Or upload a file") # uploaded_files jest katalogiem tworzonym wzgledem MEDIA_ROOT
  email = models.EmailField(max_length = 20,blank= True,help_text="Send visualization to (optional):") #blank=True, aby to pole moglo pozostac puste przy wypelnianiu pol
  unixtime = models.IntegerField(max_length = 40)
  choice = models.CharField(max_length=4, choices=gchoices,default = "Slow")
  linewidth = models.IntegerField(max_length = 300,default=30,blank=True,null = True,help_text="Number of aminoacids in one row in graph (default=30)")
    
class PageForm(forms.ModelForm):
  choice = forms.CharField(max_length=4, widget=forms.Select(choices=gchoices),help_text="Predict secondary structure (PSIPRED) with PSI-BLAST or without")
  
  def check_correctness(self,inputsequences,fields_name):
    try:
      from Bio import AlignIO
      from Bio import SeqIO
      import StringIO
      x = StringIO.StringIO(inputsequences)
      s = SeqIO.parse(x, "fasta")
    except:
      msg = "Input is not in FASTA format!"
      self._errors[fields_name] = self.error_class([msg])
      del cleaned_data[fields_name]
      return
    try:
      seqDict = SeqIO.to_dict(s)
    except ValueError:
      x.close()
      msg = "Input contains repeated names of sequences!"
      self._errors[fields_name] = self.error_class([msg])
      return
    if len(seqDict)==0:
      x.close()
      msg = "Input is not in FASTA format or the specified file is empty!"
      self._errors[fields_name] = self.error_class([msg])
      return
    l=0  
    for i in xrange(len(seqDict.values())):
      seqlength = len(seqDict.values()[i])
      if l==0:		l=seqlength
      else:
	if l!=seqlength:
	  x.close()
	  msg = "Sequence number %d isn't of the same length like the other" % i
	  self._errors[fields_name] = self.error_class([msg])
	  return

    if len(seqDict) < 2:
      x.close()
      msg = "Input contains less than 2 sequences!"
      self._errors[fields_name] = self.error_class([msg])
      return
    n=0
    for i in xrange(len(seqDict.values()[0].seq)):
      for key in seqDict.keys():
	if seqDict[key][i]!="-":
	  n+=1
    if n==0:
      x.close()
      msg = "Wrong MSA. There is no non-gap letter on position %d!" % i
      self._errors[fields_name] = self.error_class([msg])
      return

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

    if sequences!="" and upload_file==None:
      self.check_correctness(sequences,"sequences")
     # del cleaned_data["sequences"]
	
    if upload_file!=None and sequences=="":
      content = upload_file.read()
      self.check_correctness(content,"upload_file")
      #del cleaned_data["upload_file"]
      
    if linewidth < 30:
      msg = 'Linewidth must be < 30!'
      self._errors["linewidth"] = self.error_class([msg])
      del cleaned_data["linewidth"]
      
      
       
    return self.cleaned_data

  class Meta:
    model = Page
    fields = ('sequences',"upload_file", 'email',"linewidth","choice") # kolejnosc wyswietlania na stronie
   
