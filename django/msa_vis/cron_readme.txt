W terminalu należy wpisać "crontab -e" i wpisać w edytorze "0 2 * * * python /tu_sciezka_do_katalogu/msa_vis/remove.py", zapisać zmiany i wyjść. Cron codziennie o 2 w nocy uruchomi skrypt remove.py, który usuwa obrazki oraz rekordy w bazie starsze niż np. tydzień (w kodzie są też zakomentowane inne opcje - dzień, 2 tyg, miesiąc).
Można podać inną datę uruchamiania skryptu zmieniając pierwsze 5 argumentów, które oznaczają kolejno: minuty, godziny, dni miesiąca, miesiące, dni tygodnia

W razie kłopotów polecam zapisanie outputu (razem z errorem) do pliku. W edytorze crontab trzeba wpisać: 
"0 2 * * * python /tu_sciezka_do_katalogu/msa_vis/remove.py >> cron_error.txt 2>&1"
Uwaga: do pliku trzeba podać ścieżkę bezwzględną (w powyższej wersji cron_error.txt zostanie zapisany w katalogu domowym)
