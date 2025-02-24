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
import tqdm
from astropy.io import fits
from astropy.wcs import WCS

parser = argparse.ArgumentParser(description="Create line and physical parameter maps from ALFA or ALFA+NEAT analysis of data cubes. Results are written to [PREFIX]linemap.fits and [PREFIX]resultmap.fits.")
parser.add_argument("--original", default=None, required=False, help="FITS file originally analysed", nargs=1, dest="original")
parser.add_argument("directory", nargs="?", default="./", help="Directory containing output files from ALFA/NEAT")
parser.add_argument("--prefix", default="", required=False, help="prefix for output file names")
parser.add_argument("--lines-only",action='store_true',required=False,help="Only make line maps, not result maps")
parser.add_argument("--results-only",action='store_true',required=False,help="Only make result maps, not line maps")

args = parser.parse_args()

# get files

pixelfiles=glob.glob(args.directory+"/*pixel_*_*_fit.fits")
if len(pixelfiles)==0:
    print("no files found in "+args.directory)
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

if not os.path.exists(original):
  print("original file "+original+" not found. specify with --original")
  sys.exit()

print("taking WCS information from %s..."%original)

# search for a valid WCS

validwcs=False
origdata=fits.open(original)
for ext in origdata:
  if "CRPIX1" in ext.header:
    wcs = WCS(ext.header).celestial
    hdr = wcs.to_header()
    print("found WCS in extension %s"%ext.name)
    validwcs=True
    break

# fail if not found

if not validwcs:
  print("couldn't find a valid WCS")
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

# apply mapping options

if args.lines_only and args.results_only:
  print("--lines-only and --results-only are mutually exclusive (obvs)")
  sys.exit()
else:
  maplines=maplines and not args.results_only
  mapresults=mapresults and not args.lines_only

# get the fluxes

infotext=[]
if maplines:
  infotext.append(str(nlines)+" lines")
if mapresults:
  infotext.append(str(nresults)+" quantities")

print("mapping "+" and ".join(infotext)+"...")

i=0
for pixelfile in tqdm.tqdm(pixelfiles):

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
#                print(pixelfile)
                continue

    if mapresults:
        resultstable=hdus["RESULTS"].data
        for j in range(nresults):
            resultsmap[j][x][y]=resultstable[j][1]

# remove all negative values from lines image

linesmap=numpy.where(linesmap<0,0,linesmap)

# write the files
hdulist=fits.HDUList()
primaryhdu=fits.PrimaryHDU(header=hdr)
primaryhdu.header["history"] = " ".join(sys.argv)
hdulist.append(primaryhdu)

if maplines:
    print("writing file "+args.prefix+"linemap.fits...")
    for i in range(len(linelist)):
        print("%s: %i pixels mapped out of %i"%(linelist[i],numpy.count_nonzero(linesmap[i]),len(pixelfiles)))
#        if numpy.count_nonzero(linesmap[i])>0.8*len(pixelfiles):
        hdu=fits.ImageHDU(linesmap[i][:][:],header=hdr,name=str(linelist[i]))
        hdulist.append(hdu)
#        else:
#          print("didn't write extension for %s, more than 80%% non-zero"%linelist[i])
    hdulist.writeto(args.prefix+"linemap.fits",overwrite=True)

if mapresults:
    print("writing file "+args.prefix+"resultmap.fits...")
    for i in range(len(resultlist)):
        hdu=fits.ImageHDU(resultsmap[i][:][:],header=hdr,name=str(resultlist[i]))
        hdulist.append(hdu)
    hdulist.writeto(args.prefix+"resultmap.fits",overwrite=True)
