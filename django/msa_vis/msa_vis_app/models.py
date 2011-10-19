from django.db import models
from django.forms import ModelForm
from django import forms

gformat = (('fasta', 'FASTA'),
	   ('clustal','ALN/ClustalW2'),
	   ('stockholm','Pfam/Stockholm'),
	   ('phylip','Phylip'),
	   )

gchoices= (
	('Slow', 'Slower - Run PSI-BLAST'),
        ('Fast', 'Faster - without PSI-BLAST'),
        )
  
class Page(models.Model):
  sequences = models.TextField(max_length = 400,blank= True,help_text="Enter a sequence alignment (max sequences length = 400)")
  upload_file = models.FileField(upload_to = "uploaded_files", blank= True, help_text="Or upload a file")
  email = models.EmailField(max_length = 60,blank= True,help_text="Send visualization to (e-mail address; optional)")
  unixtime = models.IntegerField(max_length=50)
  choice = models.CharField(max_length=4, choices=gchoices,default = "Slow")
  linewidth = models.IntegerField(max_length = 1000,default=100,blank=True,null = True,help_text="Number of amino acids in one row in graph")
  format = models.CharField(max_length=20, choices=gchoices,default = "FASTA")
    
class PageForm(forms.ModelForm):
  choice = forms.CharField(max_length=4, widget=forms.Select(choices=gchoices),help_text="Predict secondary structure (PSIPRED) with PSI-BLAST or without")
  format = forms.CharField(widget=forms.Select(choices=gformat),help_text="Select input sequence format")
  
  
  def check_correctness(self,inputsequences,fields_name,format):
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
    except ValueError:
      x.close()
      msg = "Input contains repeated names of sequences"
      self._errors[fields_name] = self.error_class([msg])
      return
    except:
      x.close()
      msg = "Selected wrong input format!"
      self._errors[fields_name] = self.error_class([msg])
      return
    if len(seqDict)==0:
      x.close()
      msg = "Input is not in %s format or the specified file is empty!" % (format)
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
      self.check_correctness(sequences,"sequences",format)
      #del cleaned_data["sequences"]
	
    if upload_file!=None and sequences=="":
      content = upload_file.read()
      self.check_correctness(content,"upload_file",format)
      #del cleaned_data["upload_file"]
      
    if linewidth < 1:
      msg = 'Linewidth must be >= 1!'
      self._errors["linewidth"] = self.error_class([msg])
      del cleaned_data["linewidth"]
      
    
    return self.cleaned_data

  class Meta:
    model = Page
    fields = ('format','sequences',"upload_file", 'email',"linewidth","choice")
   