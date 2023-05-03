# GLM_bolide_tools

## dependencies
```
numpy pandas astropy wget
```

## download
Use `python download_and_process_GLM_bolides.py`.

This will create a directory structure (`./glm_bolides_events` by default), save the general catalogue to both json and votable, and download all light curves / navigated positions as CSV files.

## convert
`python GLM_to_GFE.py [PATH/TO/FILE.csv]`

This will convert a CSV file to the Global Fireball Exchange astrometric standard (https://github.com/UKFAll/standard).

It tries to de-navigate the lat-long values assuming the lightning ellipsoid.
Tested on one event, result seem to integrate well when solving with other ground-based astrometry.
But needs more validation.

