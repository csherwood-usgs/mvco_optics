# mvco_optics
Routines, notes, and results for review paper on particle optics


`plot_pfa_ts`
* Loads the time series data `ustar_av.mat` created by `plot_ustar.m`
* Loads LISST size bins from `lisst_av`
* Loads a bunch of `*_pfa.mat` files
* Interpolates the results onto the u* time base


`_pfa.mat' files are time series of profile fit parameters
