#!/usr/bin/python

# this is a script to turn ALFA and NEAT FITS output into image maps
# rw@nebulousresearch.org
# License: GPL v3
# todo: add --no-wcs option in case original file not available
#       handle case when section of original file was analysed

import argparse
import glob
import math
import numpy
import os
import sys
from astropy.io import fits
from astropy.wcs import WCS

parser = argparse.ArgumentParser(description="Create line and physical parameter maps from ALFA or ALFA+NEAT analysis of data cubes")
parser.add_argument("--original", default=None, required=False, help="FITS file originally analysed", nargs=1, dest="original")
args = parser.parse_args()

# get files

pixelfiles=glob.glob("./*_*_*_fit.fits")
if len(pixelfiles)==0:
    print("no files")
    sys.exit()
else:
    print("mapping "+str(len(pixelfiles))+" pixels...")

header=fits.getheader(pixelfiles[0])

# original file from command line if specified, otherwise try to get it from the header

if args.original:
  original=args.original[0]

else:
  original=""
  history="".join(header["history"]).split(" ")

  for s in history:
    if ".fits" in s:
      original=s
      break

  if original=="":
    print("couldn't deduce original file from header information. specify with --original")
    sys.exit()

# if the original file exists, get the WCS

try:
  hdu = fits.open(original)[1].header
  wcs = WCS(hdu).celestial # Import the WCS header
  hdr = wcs.to_header()
except:
  print("original file "+original+" not found. specify with --original")
  sys.exit()

# get the dimensions. filenames have y coordinate first

xpix=[]
ypix=[]

for pixelfile in pixelfiles:
    xpix.append(int(pixelfile.split("_")[-2])-1)
    ypix.append(int(pixelfile.split("_")[-3])-1)

# work out what to map

hdu=fits.open(pixelfiles[0])
maplines=True
mapresults=True

try:
    table=hdu["LINES"].data
    nlines=len(table)
    linesmap=numpy.zeros(shape=(nlines,max(xpix)+1,max(ypix)+1))
    wlenlist=table["WlenRest"]
    ionlist=table["Ion"]
    linelist=[ionlist[x]+" "+str(wlenlist[x]) for x in range(len(ionlist))]
except:
    print("no lines")
    nlines=0
    maplines=False

try:
    table=hdu["RESULTS"].data
    nresults=len(table)
    resultsmap=numpy.zeros(shape=(nresults,max(xpix)+1,max(ypix)+1))
    resultlist=table["Quantity"]
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
    hdulist=fits.HDUList()
    for i in range(len(linelist)):
        if numpy.any(linesmap[i]>0):
          hdu=fits.ImageHDU(linesmap[i][:][:],header=hdr,name=str(linelist[i]))
          hdulist.append(hdu)
    hdulist.writeto("linemap.fits",overwrite=True)

if mapresults:
    hdulist=fits.HDUList()
    for i in range(len(resultlist)):
        hdu=fits.ImageHDU(resultsmap[i][:][:],header=hdr,name=str(resultlist[i]))
        hdulist.append(hdu)
    hdulist.writeto("resultmap.fits",overwrite=True)
