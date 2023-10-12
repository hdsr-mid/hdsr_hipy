Deze environment-directory bevat de volgende sub-dirs en bestanden:
- `HYDROMEDAH_python38`: de volledige python-environment voor hydsrhipy
- `wheels`: de packages uit `environment_export.yml` gedownload als wheels
- `environment.yml`: de Python-environment waarmee `HYDROMEDAH_python38` is aangemaakt met micromamba op 11-10-2023
- `environment_export.yml`: de environment-specificatie zoals geexporteerd uit `HYDROMEDAH_python38`
- `environment_docs.yml`: de environment voor het publiceren van documentatie

De exacte `HYDROMEDAH_python38` environment is na te maken met [micromamba](https://prefix.dev/docs/mamba/overview#installation):
```
micromamba env create -f environment_export.yml
```

Niet aanbevolen, maar als alternatief kan de installatie worden uitgevoerd door alle wheels te installeren met pip:
```
pip install <filename>.whl --no-deps
```

Een beter idee is om de `HYDROMEDAH_python38` folder in een bestaande Mamba/Conda te zetten.
