# mvco_optics
Routines, notes, and results for review paper on particle optics

### Data files
`advpr_mvco11_05-Jan-2012.mat` - Processed ADV pair data, with velocities, wave parameters, and stresses
`suspsed_ba20_rstrim_crs.mat` - Profile data structure, with fields like this:
```
ba = 
  struct with fields:

         ADVagc: [12×2441 double]
            OBS: [12×2441 double]
       tranattn: [12×2441 double]
              u: [12×2441 double]
              v: [12×2441 double]
          ADVdn: [12×2441 double]
          OBSdn: [12×2441 double]
         trandn: [12×2441 double]
           ADVz: [12×2441 double]
           OBSz: [12×2441 double]
          tranz: [12×2441 double]
          ADVbe: [0.0500 0.2500 0.4500 0.6500 0.8500 1.0500 1.2500 1.4500 1.6500 1.8500 2.0500 2.2500 2.4500]
      ADVfields: {13×1 cell}
      LISSTattn: [12×2442 double]
      LISSTD16v: [12×2442 double]
      LISSTD50v: [12×2442 double]
      LISSTD84v: [12×2442 double]
       LISSTsvv: [12×2442 double]
     LISSTSzmnv: [12×2442 double]
      LISSTvtot: [12×2442 double]
      LISSTD16a: [12×2442 double]
      LISSTD50a: [12×2442 double]
      LISSTD84a: [12×2442 double]
       LISSTsva: [12×2442 double]
     LISSTSzmna: [12×2442 double]
      LISSTatot: [12×2442 double]
    LISSTfinesv: [12×2442 double]
    LISSTmicrov: [12×2442 double]
    LISSTmacrov: [12×2442 double]
    LISSTfinesa: [12×2442 double]
    LISSTmicroa: [12×2442 double]
    LISSTmacroa: [12×2442 double]
     LISSTvconc: [12×2442×32 double]
     LISSTaconc: [12×2442×32 double]
        LISSTsz: [32×1 double]
        LISSTdn: [12×2442 double]
         LISSTz: [12×2442 double]
        LISSTbe: [0.0500 0.2500 0.4500 0.6500 0.8500 1.0500 1.2500 1.4500 1.6500 1.8500 2.0500 2.2500 2.4500]
    LISSTfields: {26×1 cell}
    ```
    

`_pfa.mat` files are time series of profile fit parameters

`plot_pfa_ts.m` - prepares data for `fdyn.m`
* Loads the time series data `ustar_av.mat` created by `plot_ustar.m`
* Loads LISST size bins from `lisst_av`
* Loads a bunch of `*_pfa.mat` files
* Interpolates the results onto the u* time base, then onto two-hour time base

`fdyn.m` - various time series plots of floc parameters

