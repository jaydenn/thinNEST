# thinNEST
A lightweight and flexible code for running NEST simulations. thinNEST is modelled off the testNEST code, but with a few added features: 
- simulate arbitrary spectra supplied via an input file
- optional inclusion of the migdal effect
- optional calculation of the s2-width (requires modified nest source)
- changable detector and analysis parameters

With the goal of running fast and without having to re-compile to change options.

## Dependencies

- NEST: https://github.com/NESTCollaboration/nest
- GSL (best installed via a package manager)

## Compilation

First edit the Makefile to set the location of your NEST library directory (by default it is assumed to be in the parant dir).
Then type 'make'.

## Running

Usage: thinNEST [OPTIONS]  

Options:
        -v, --verbose   verbose output
        -q, --outputQuanta      controls what is saved in output
        -t, --timing    include s2-width calculation
        -m, --migdal    include the migdal effect for NR
        -M, --migdal    include the migdal effect for NR (optimized, relative NR/Mig rate will be unphysical)
        -n, --numEvents N       simulate N events
        -a, --analysis analysis_file
        -d, --detector detector_file
        -s, --spectrum {NR,ER,file=spectrum_file}
        -e, --eMin X    starts the spectrum at X keV, defaults to -1 (min in spectrum file)
        -E, --eMax X    ends the spectrum at X keV, defaults to -1 (max in spectrum file)
        -o, --output output_file, defaults to stdout
        -f, --field V   set E-field to V volts/cm
