import sys
import os
import csv

from argparse import ArgumentParser
from collections import defaultdict

from pprint import pprint

import numpy as np

# ====================================================================


def main():
    
    # Parse arguments from command line
    parser = ArgumentParser()
    parser.add_argument("--plusfiles",nargs="+",
                        help="Tabular file with read 5' ends for the + strand. If multiple files are provided, read counts will be pooled.")
    parser.add_argument("--minusfiles",nargs="+",
                        help="(Optional) Tabular file with read 5' ends for the - strand. If multiple files are provided, read counts will be pooled.")
    parser.add_argument("--left",type=int,
                        help="left-most position of the target window")
    parser.add_argument("--right",type=int,
                        help="right-most position of the target window")
    parser.add_argument("--datadir",default="",
                        help="Where the data files are stored, if not working directory")
    parser.add_argument("-o","--outfile",
                        help="Output file for peaks")

    args = parser.parse_args()

    # Make sure the directory name ends with a '/'
    if args.datadir and args.datadir[-1] != "/":
        args.datadir = args.datadir + "/"
    if not args.datadir:
        args.datadir=""

    # Call peaks on the reads mapping to the (+) strand
    peaks =  call_peaks([ args.datadir+fname for fname in args.plusfiles ],
                        args.left,args.right,"+")
    if args.minusfiles:
        peaks.extend( call_peaks([ args.datadir+fname for fname in args.minusfiles ],
                                 args.left,args.right,"-") )

    # Write peaks to output file
    with open(args.outfile,'w') as out:
        header = [ "chrom","pos","strand","Z","std","reads"]
        out.write( "\t".join(header) )

        # Sort the entries by sequence and position
        peaks = sorted(peaks,
                       key=lambda x: (x[0],x[1]) )
        for result in peaks:
            result = [ str(k) for k in result ]
            out.write( "\n" + "\t".join(result) )
                
        
# ====================================================================

def call_peaks(filename,
               left,right,
               strand):

    results = [ ] # List of results,
                  # entries will be (chrom,pos,strand,Z,Zplus,Zminus)

    # Load reads
    pooled = len(filename) > 1
    reads  = load_reads(filename,pooled)

    # Loop over target sequences
    for chrom,chromreads in reads.items():

        # Select the target window, calculate total reads in window, then
        # normalize reads at each position
        chrom_positions = [ i for i in range(left,right+1) ]
        total = float(sum( [chromreads[p] for p in chrom_positions] ))
        myreads = { p : chromreads[p]/total for p in chrom_positions }

        # Skip regions with extremely low coverage
        if total < 10:
            continue

        # For each position in the window
        for pos in chrom_positions:

            # Fetch the values at the position of interest, and the
            # background set
            stop = myreads[pos]
            bg = [ myreads[i] for i in chrom_positions
                   if i not in [pos,pos-1,pos+1] ]

            # Calculate Z score and store the result
            Z,std = calculate_Z(stop,bg)
            result = ( chrom,pos,strand,          # Three fields for the position
                       np.round(Z,decimals=3),    # The Z-score
                       np.round(std,decimals=3),  # The std component of the Z-score
                       total )                    # and the window coverage
            results.append( result )

    return results

# ====================================================================


def load_reads(infiles,pooled=False):

    reads = defaultdict( lambda:defaultdict(float) )

    for fname in infiles:

        f = open(fname)

        # Read ends from .wig
        fileinfo = fname.split(".")
        if fname.split(".")[-1] == "wig":
            f.readline()
            while True:
                # Read a line, make sure it's not EOF
                line = f.readline().strip().split()
                if not line:
                    break
                # Parse the chromosome line when I find it
                if line[0] == "variableStep":
                    chrom = line[1].split("=")[-1]
                # If the line contains position info, store it
                # for the chromosome I saw most recently
                else:
                    pos,count = line
                    reads[chrom][int(pos)] += float(count)

        # Read ends from .txt file
        else:
            for r in csv.DictReader(f,delimiter="\t",
                                    fieldnames=["chrom","pos","count"]):
                count = float(r['count'])
                if count > 0:
                    reads[ r['chrom'] ][ int(r['pos']) ] += count
        f.close()

    return reads

        
def calculate_Z(stop,bg):
    mu = np.mean(bg)
    sigma = np.std(bg)
    Z = (stop - mu) / sigma
    return Z,sigma


# ====================================================================

if __name__=="__main__":
    main()
