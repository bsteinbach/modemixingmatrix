#!/usr/bin/env python

import argparse
import os
import healpy as hp
import numpy as np

def generate_powerspectrum(mask,lmax):
	print "Reading mask file"
	m = hp.read_map(mask)
	clfn = mask + '.cl'
	nside = hp.npix2nside(m.size)
	print "Generating power spectrum"
	cl = hp.anafast(m,lmax=lmax)
	np.savetxt(clfn,cl)
	return clfn

def run_modemixing(outfn,clfn,lmax):
	print "Calculating mode mixing matrix"
	cmd = './modemix %s %s %d'%(clfn,outfn,lmax)
	print cmd
	os.system(cmd)

def main():
	parser = argparse.ArgumentParser('Compute mode mixing matrix for a healpix mask')
	parser.add_argument('--mask',required=True,help='Name of healpix mask file')
	parser.add_argument('-o',help='Name of output mode mixing matrix',required=True)
	parser.add_argument('-l','--lmax',help='lmax of output mode mixing matrix',type=int,required=True)
	args = parser.parse_args()

	clfn = generate_powerspectrum(args.mask,args.lmax)
	run_modemixing(args.o,clfn,args.lmax)

if __name__=='__main__':
	main()
