Instalacja django-extensions.

W Ubuntu wystarczy zainstalowac pakiet python-django-extensions. Zeby dzialy niektore fragmenty django-extensions potrzebne sa jeszcze inne pakiety, np. zeby dzialalo debugowanie przez strone www trzeba doinstalowac python-werkzeug.

W settings.py nalezy dodac 'django_extensions' do INSTALLED_APPS. Wtedy w manage.py beda dostepne dodatkowe komendy, np. ./manage.py runserver_plus (uruchamia serwer z debugowaniem przez werkzeug).