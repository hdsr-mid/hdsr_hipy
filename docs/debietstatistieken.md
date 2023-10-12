# Debietstatistieken ten behoeve van thermische energie uit oppervlaktewater
HDSR Hydrologische informatieproducten - debietstatistiek

Ten behoeve van het project 'Debietstatistieken ten behoeve van thermische energie uit oppervlaktewater' is HDSRHIPY uitgebreid met de volgende features:

- Conversie van laterale afvoer naar een SOBEK LATERAL.DAT-bestand
- Functie die statistieken berekend voor een op te geven periode in het jaar op basis van Sobek resultaten




<br><br>

## Conversie van laterale afvoer naar een SOBEK LATERAL.DAT-bestand
Uit de HYDROMEDAH run komt een CSV-bestand met de afwateringseenheden als kolommen en de tijdtappen als rijen. Deze informatie gebruiken we om een LATERAL.DAT te vullen voor Sobek.
Hiervoor is de volgende informatie nodig:
- een 'leeg' LATERAL.DAT bestand zonder tijdreeksen, maar met de juiste ID's van de lateralen (meegeleverd in de 'templates'-map);
- GIS-bestanden die zijn geexporteerd uit SOBEK (3b_nod.shp en 3b_link.shp);
- een koppeltabel tussen Hydromedah-SWNR-units en SOBEK-IDS van afwateringseenheden,.

Vervolgens wordt voor elke lateral in Sobek de juiste kolom opgezocht in de CSV en wordt de reeks in het juiste format geschreven.

Het kan voorkomen dat numerieke problemen optreden. Het is daarom mogelijk een dictionary mee te geven met een koppeling tussen reach-segmenten en lateralen. Deze zijn in GIS viseel gekoppeld en via de volgende stappen tot stand gekomen:
- uitlezen FLOWANAL.HIS en bepalen op welke reach-segmenten de meeste iteraties veroorzaken;
- in GIS uitzoeken welke lateralen daar net boventrooms van of op liggen. 
Op alle lateralen die voorkomen in de dictionary wordt een anti-droogvaldebiet van 0.01 m3/s toegepast. De LATERAL.DAT is vervolgens gebruikt om de periode in SOBEK door te rekenen. Dit is niet opgenomen in HDSRHIPY omdat het te tijdrovend is.

Als voorbeeld: de dict "{}'H013921_1':['PG0599-1','PG0034-1','PG2177-1']}" geeft aan dat aan reach 'H013921_1' drie lateralen gekoppeld zijn: PG0599-1, PG0034-1, en PG2177-1.

<br><br>

## Functie die statistieken bereket voor een op te geven periode in het jaar op basis van Sobek resultaten

Deze functie is inbegrepen in hdsrhipy en is een verkorte weergave van het eerder opgeleverde notebook 'debietstatistieken.ipynb'. Dat notebook is ook nu weer inbegrepen.

In hdsrhipy moet eerst een 'flowstats' object worden aangemaakt. 

```{python}
flowstats = FlowStats(sobek_path=hisfiles, template_path=template_path, export_path=os.path.join(export_path, 'flowstats.shp'))
```

Hierbij is:
- sobek_path: het pad naar de Sobek uitvoer, tot en met de .lit-map;
- template_path: de plek waar de shape met reachsegmenten, geexporteerd uit Sobek, staat;
- export path: de plek waar de resulterende shapefile wordt geschreven.

Vervolgens worden de statistieken uitgerekend met de functie "get_flowstats()": 
```{python}
flowstats.get_flowstats(cases2join=cases, period=period, overschrijdingsduur_van_debiet=overschrijdingsduur_van_debiet, relatieve_overschrijdingsduur_van_debiet=relatieve_overschrijdingsduur_van_debiet,debiet_bij_tijdpercentiel=debiet_bij_tijdpercentiel)
```

waarbij "cases" de lijst van SObek-cases is waarvan de HIS-uitvoer moet worden samengevoegd. De andere argumenten geven de periode waarover statistieken worden uitgerekend en de parameters ervan:

```{python}
# maanden in het jaar om statistiek over te berekenen. De maand- en dag van deze periode wordt gebruikt om te bepalen of de dag moet worden meegenomen
period = pd.date_range(start=pd.Timestamp('1900-05-15'),end=pd.Timestamp('1900-09-15'), freq='D')

# de 'name' wordt gebruikt als kolom in de shapefile, de waarde is het debiet (m3/s) om te bepalen hoe vaak het wordt overschreden. Het resultaat is het aantal dagen.
overschrijdingsduur_van_debiet = {'name': 'Ex5-9_005',  'value': 0.05}
# Idem, maar nu is het resultaat de fractie van de tijd dat de waarde wordt overschreden
relatieve_overschrijdingsduur_van_debiet = {'name': 'REx5-9_005',  'value': 0.05}
# en dit geeft het debiet dat x % van de tijd wordt overschreden
debiet_bij_tijdpercentiel = {'name': 'Q5-9_80%',  'value': 0.8}
```