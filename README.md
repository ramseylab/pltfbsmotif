# pltfbsmotif

**Piecewise-linear machine-learning model for motif-based prediction of transcription factor binding sites**

This repository contains the MATLAB code for a piecewise-linear, data-fusion model for
motif-based prediction of transcription factor binding sites. The code was used in
the analysis published in the 2010 article *Genome-wide histone acetylation data improve prediction of mammalian transcription factor binding sites* 
[(Ramsey et al., Bioinform., 2010)](https://doi.org/10.1093/bioinformatics/btq405).

The MATLAB code in this project was written by Stephen Ramsey [(@saramsey)](https://github.com/saramsey) except where otherwise indicated in the header comments
and except for the following code modules which were either contributed or incorporated
from public-domain sources. 

- PSWM_MotifLocator.m:  Harri Lahdesmaki (personal communication)
- readbar.m/writebar.m:  Matti Nykter (personal communication)
- uuidgen.m: Dahua Lin (from MATLAB Central)

The aforementioned authors of these contributed and public-domain software modules
are gratefully acknowledged.

This software requires the MATLAB software Stable Noisy Optimization by Branch and FIT (SNOBFIT)
by Waltraud Huyer and Arnold Neumaier and which is freely available for download from Universitat Wien
(download SNOBFIT)[http://www.mat.univie.ac.at/~neum/software/snobfit/].

This software also requires the MATLAB Statistics Toolbox.

