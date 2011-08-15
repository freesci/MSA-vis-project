from django import forms

class ContactForm(forms.Form):
    email = forms.EmailField(max_length = "20")
    #sequences = forms.TextField(max_length = "150") # help_text nie jest najladniejszym rozwiazaniem w widoku formularzu...
    #upload_file = forms.FileField(upload_to = "uploaded_files") #storage=fs # uploaded_files jest katalogiem tworzanym wzgledem MEDIA_ROOT
    #seqID = models.CharField(max_length = "40") # podobe pole do tego ze strony aln2plot - tam jesli uzytkownik nic w to pole nie wipisze, program sam wygeneruje losowe id