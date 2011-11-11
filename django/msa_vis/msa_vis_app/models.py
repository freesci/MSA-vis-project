from django.db import models
from django.forms import ModelForm
from django import forms

gformat = (('fasta', 'FASTA'),
	   ('clustal','CLUSTAL(<=2.0)'),
	   ('phylip','Phylip4|Phylip(interleaved)'),
	   )

gchoices= (
	('Slow', 'Slower - Run PSI-BLAST'),
        ('Fast', 'Faster - without PSI-BLAST'),
        )

class Page(models.Model):
  sequences = models.TextField(max_length = 400,blank= True,help_text="Enter a sequence alignment (max sequences length = 400)")
  upload_file = models.FileField(upload_to = "uploaded_files/input_files", blank= True, help_text="Or upload a file")
  email = models.EmailField(max_length = 60,blank= True,help_text="Send visualization to (e-mail address; optional)")
  unixtime = models.IntegerField(max_length=50)
  choice = models.CharField(max_length=4, choices=gchoices,default = "Slow")
  linewidth = models.IntegerField(max_length = 1000,default=100,blank=True,null = True,help_text="Number of amino acids in one row in graph")
  format = models.CharField(max_length=20, choices=gchoices,default = "FASTA")

class PageForm(forms.ModelForm):
  choice = forms.CharField(max_length=4, widget=forms.Select(choices=gchoices),help_text="Predict secondary structure (PSIPRED) with PSI-BLAST or without")
  format = forms.CharField(widget=forms.Select(choices=gformat),help_text="Select input sequence format")


  def check_correctness(self,inputsequences,fields_name,format,cleaned_data):
    try:
      from Bio import AlignIO,SeqIO
      import StringIO
      x = StringIO.StringIO(inputsequences)
      s = SeqIO.parse(x, format)
    except:
      msg = "Input is not in %s format!" % (format)
      self._errors[fields_name] = self.error_class([msg])
      del cleaned_data[fields_name]
      return
    try:
      seqDict = SeqIO.to_dict(s)
      x.close()
    except ValueError:
      msg = "Input contains repeated names of sequences"
      self._errors[fields_name] = self.error_class([msg])
      x.close()
      return
    except:
      msg = "Selected wrong input format!"
      self._errors[fields_name] = self.error_class([msg])
      return
    if len(seqDict)==0:
      msg = "Input is not in %s format or the specified file is empty!" % (format)
      self._errors[fields_name] = self.error_class([msg])
      return
    l=0
    for i in xrange(len(seqDict.values())):
      seqlength = len(seqDict.values()[i])
      if l==0:		l=seqlength
      else:
	if l!=seqlength:
	  msg = "Sequence number %d isn't of the same length like the other" % (i+1)
	  self._errors[fields_name] = self.error_class([msg])
	  return
    if len(seqDict) < 2:
      msg = "Input contains less than 2 sequences!"
      self._errors[fields_name] = self.error_class([msg])
      return
    n=0
    for i in xrange(seqlength):
      for key in seqDict.keys():
	if seqDict[key][i]!="-" and seqDict[key][i]!="X" and seqDict[key][i]!="x":
	  n+=1
    if n==0:
      msg = "Wrong MSA. There is no non-gap letter on position %d!" % i
      self._errors[fields_name] = self.error_class([msg])
      return

    cleaned_data["seqlength"]=seqlength



  def clean(self):
    cleaned_data = self.cleaned_data

    sequences = cleaned_data["sequences"]
    linewidth = cleaned_data["linewidth"]
    upload_file = cleaned_data["upload_file"]
    format = cleaned_data["format"]

    if sequences=="" and upload_file==None:
      msg = 'You must specify an input source!'
      self._errors["sequences"] = self.error_class([msg])
      del cleaned_data["sequences"]

    if sequences!="" and upload_file!=None:
      msg = 'Select alignment or upload a local file'
      self._errors["sequences"] = self.error_class([msg])
      del cleaned_data["sequences"]

    if sequences!="" and upload_file==None:
      self.check_correctness(sequences,"sequences",format,cleaned_data)
      #del cleaned_data["sequences"]

    if upload_file!=None and sequences=="":
      content = upload_file.read()
      self.check_correctness(content,"upload_file",format,cleaned_data)
      #del cleaned_data["upload_file"]

    if linewidth < 1:
      msg = 'Linewidth must be >= 1!'
      self._errors["linewidth"] = self.error_class([msg])
      del cleaned_data["linewidth"]

    return self.cleaned_data

  class Meta:
    model = Page
    fields = ('format','sequences',"upload_file", 'email',"linewidth","choice")
