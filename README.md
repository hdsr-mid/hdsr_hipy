# hdsrhipy
HDSR Hydrologische informatieproducten


## Installatie Python package

1.	Installeer een Anaconda of Miniconda Python distributie:

•	https://www.anaconda.com/products/individual
•	https://docs.conda.io/en/latest/miniconda.html


2.	Ga naar de URL: 

https://github.com/RuudHurkmans/hdsrhipy/blob/master/environment.yml

En kopieer de inhoud naar een lokaal bestand, eveneens genaamd ‘environment.yml’. Op dit moment is de ‘hdsrhipy’- respository nog ‘private’ en is dit bestand nog niet toegankelijk zonder uitnodiging. Projectteamleden zijn uitgenodigd. Als bovenstaande URL niet werkt, stuur dan een e-mail naar hurkmans@hkv.nl.

3.	Open een Ananconda-prompt door op start te klikken en ‘anaconda in te typen’.

4.	Navigeer in de prompt naar de locatie waar het in stap 2 gekopieerde bestand staat:

cd <pad naar environment.yml>

5.	Type het commando “conda env create -f environment.yml” en druk op <ENTER>.

Nu wordt een Python-omgeving genaamd ‘hdsrhipy’ aangemaakt met daarin alle benodigde afhankelijkheden. Het pakket hdsrhipy zelf hoort daar nog niet bij. In een later stadium kan dit eventueel worden opgewaardeerd tot een openbaar en installeerbare Pyton-bibliotheek, maar voor nu dient de code handmatig te worden gedownload of gekloond vanaf Github. 

6.	Ga naar de URL:

https://github.com/RuudHurkmans/hdsrhipy

En download de code als zip-bestand via “Download ZIP”. Pak de code op de gewenste locatie uit. Merk op dat op dit moment hdsrhipy op dit moment nog niet openbaar is en niet toegangkelijk zonder Github account. Hiet ZIP-bestand is nu bijgevoegd bij de oplevering. Pak deze dus uit op de gewenste locatie. De laatste stap is het registren van bovenstaand environment (hdsrhipy) aan Jupyter notebook. Voer daartoe dit commando uit:

python -m ipykernel install --user --name=hdsrhipy

Nu zijn alle benodigde onderdelen geïnstalleerd en kan een Jupyter-notebook worden opgestart om de code uit te voeren.

4.2	Gebruik Jupyter-notebook

Voer, om het notebook op te starten, de volgende commando’s uit in de Anaconda prompt.

-	conda activate hdsrhipy (de aangemaakte Python-omgeving wordt actief).
-	cd /D <PAD NAAR UITEPAKTE ZIP>\notebooks
-	jupyter-lab

Nu opent in de browser het betreffende notebook. Om de code daarin te gebruiken, dienen in de bovenste cel een paar paden te worden aangepast:
-	Package_path: het pad naar de map waarin de zip is uitgepakt;
-	Data_path: het pad waarin de Hydromedah uitvoer te vinden is;
-	Export_path: het pad waarin de figuren en de resulterende GIS-bestanden worden geplaatst. 
Met CNTR-ENTER kan een cel worden uitgevoerd. Verdere uitleg is te vinden in het notebook zelf en bijlage A.
