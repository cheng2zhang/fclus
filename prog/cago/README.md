# Flat-histogram sampling for alpha-carbon Go model

## Files

 File         | Description
--------------|------------------------------
rmsdmd.c      | molecular dynamics with hybrid MC, flat histogram along RMSD
rmsdmc.c      | Monte Carlo, flat histogram along RMSD
ncmd.c        | molecular dynamics with hybrid MC, flat histogram along the number of contacts
ncmc.c        | Monte Carlo, flat histogram along the number of contacts


## Usage

### rmsdmd

Explicit HMC, regular run

```
./rmsdmd
```

Implicit HMC, regular
```
./rmsdmd --ihmc
```
