Calculate mode mixing matrices for unbiased masked 2D power spectrum analysis on a healpix grid in a MASTER style pipeline.  This code is a direct implementation of the appendix of Hivon et al "MASTER of the CMB Anisotropy Power Spectrum: A Fast Method for Statistical Analysis of Large and Complex CMB Data Sets" http://arxiv.org/abs/astro-ph/0105302.

Dependencies:
gcc
numpy
healpy

Compilation:
make

Usage:
./modemix.py --mask maskfn -o mll.txt -l 1536

lmax is typically three times nside.
