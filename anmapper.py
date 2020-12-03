#!/usr/bin/python

# this is a script to turn ALFA and NEAT FITS output into image maps
# rw@nebulousresearch.org
# todo: copy WCS from input FITS
# License: GPL v3

import glob
import math
import numpy
import sys
from astropy.io import fits

# get files

pixelfiles=glob.glob("./*_*_*_fit.fits")
if len(pixelfiles)==0:
    print("no files")
    sys.exit()
else:
    print("mapping "+str(len(pixelfiles))+" pixels...")
# get the dimensions

xpix=[]
ypix=[]

for pixelfile in pixelfiles:
    xpix.append(int(pixelfile.split("_")[-3]))
    ypix.append(int(pixelfile.split("_")[-2]))

# work out what to map

hdu=fits.open(pixelfiles[0])
maplines=True
mapresults=True

try:
    table=hdu["LINES"].data
    nlines=len(table)
    linesmap=numpy.zeros(shape=(nlines,max(xpix)+1,max(ypix)+1))
except:
    print("no lines")
    nlines=0
    maplines=False

try:
    table=hdu["RESULTS"].data
    nresults=len(table)
    resultsmap=numpy.zeros(shape=(nresults,max(xpix)+1,max(ypix)+1))
except:
    print("no results")
    nresults=0
    mapresults=False

hdu.close()

# get the fluxes

print("mapping "+str(nlines)+" lines and "+str(nresults)+" quantities...")

i=0
for pixelfile in pixelfiles:

    x=xpix[i]
    y=ypix[i]
    i=i+1

    hdus=fits.open(pixelfile)

    if maplines:
        linestable=hdus["LINES"].data
        for j in range(nlines):
            try:
                linesmap[j][x][y]=linestable[j][2]
            except:
                print(pixelfile)
                continue

    if mapresults:
        resultstable=hdus["RESULTS"].data
        for j in range(nresults):
            resultsmap[j][x][y]=resultstable[j][1]

# remove all negative values from lines image

linesmap=numpy.where(linesmap<0,0,linesmap)

# write the files

print("writing file...")
if maplines:
    fits.writeto("linemap.fits",linesmap,overwrite=True)

if mapresults:
    fits.writeto("resultmap.fits",resultsmap,overwrite=True)
