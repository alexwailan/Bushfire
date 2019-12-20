#!/usr/bin/env python3


#Created: 27.11.18 - Alexander Wailan 



import os
import sys
import gzip
import argparse
import pandas as pd
import numpy as np
import subprocess
import pathlib
from pathlib import Path
import shutil

##########################################
# Check if dependent programs are loaded #
##########################################

dependencies = [
    'remove_blocks_from_aln.py',
    'snp-sites',
    'iqtree',
    'pyjar.py',
    'potplant.py',
    'purge.py',
    'germie.py'
]


def depend_check(dependencies):
    all_d = []
    for d in dependencies:
        if shutil.which(d, mode=os.F_OK | os.X_OK, path=None) is not None:
            print("%s has been found!" %d)
            all_d.append('TRUE')
        elif shutil.which(d, mode=os.F_OK | os.X_OK, path=None) is None:
            print("Unable to find %s." %d)
            all_d.append('FALSE')
    return all_d

##########################################
# How one communicates an error
##########################################

def ErrorOut(error):
    print ("\nError: ", error)
    print ("\nThat's pretty sad face. Double check all inputs using -h or --help. Or call me ... maybe?")
    print ()
    sys.exit()




##########################################
# Function to Get command line arguments #
##########################################

def getargv():
    usage = 'bushfire.py [options] reference ids'
    description='Run Bushfire. A program to parse a multi-fasta core SNV alignment file to a ML tree'
    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument('aln_file', help='Core SNV alignment file', metavar='FILE',type=str, nargs='?')
    parser.add_argument('masking_file', help='The masking files with regions defined in an EMBL style tab-delimited file.', metavar='FILE', type=str, nargs='?')
    parser.add_argument('tree', action="store", help='Name of output JAR tree.', metavar='PREFIX', type=str, nargs='?',default='tree')
    parser.add_argument('-d',	'--dirpath', help='Specify input directory containing files. End with a forward slash. Eg. /temp/fasta/', metavar='DIR', type=str, nargs='?', default=os.getcwd()+'/') 
    parser.add_argument('-o',	'--outdir', help='Specify output directory. End with a forward slash. Eg. /temp/fasta/; Default to use current directory.', metavar='DIR', type=str, nargs='?', default=os.getcwd()+'/')       
    return parser.parse_args()

##########################################
    # Functions for the main program #
##########################################

def purge():

    command = "purge.py %s %s"%(idir+afile,idir+mfile) #building the command to be run

    p = subprocess.call(command, shell=True,    stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    gzip_f = Path(odir+i+'/'+i+'.vcf.gz')

    ##########################################
    # Check if masked aln file exists (Purge out file)
    ##########################################

    mcore_f = Path(idir+'masked_core.snp_sites.aln')
    if not os.path.isfile(mcore_f):
        ErrorOut('Purge failed.')

    print("Purge successfull!")
    
def germie():

    command = "germie.py %s %s"%(idir+purgefile,tree)

    p = subprocess.Popen(command, shell=True,    stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output,error = p.communicate() #Read data from stdout and stderr. Wait for process to terminate.
    
    ##########################################
    # Check if JAR tree file exists (Germie out file)
    ##########################################

    jar_f = Path(idir+tree+'.joint.tre')

    if not os.path.isfile(jar_f):
        ErrorOut('germie failed.')    
    
def potplant():
    
    command = "potplant.py %s"%(idir+JARtreefile)

    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    output,error = p.communicate() #Read data from stdout and stderr. Wait for process to terminate.
    
    ##########################################
    # Check if Potplant file exists (potplant out file)
    ##########################################
    potp_f = Path(idir+tree+'.joint.tre')

    if not os.path.isfile(potp_f):
        ErrorOut('potplant failed.')


################
# Main program #
################

def main():
    #############################################################################################
    #
    #      Parse/ check the arguements
    #
    #############################################################################################

    args = getargv()

    if args.aln_file == None:
        ErrorOut('No alignment file stated.')
    elif args.masking_file == None:
        ErrorOut('No masking file stated. I need a mask!')

    idir = args.dirpath ##the working directory that holds the samples
    odir = args.outdir ##the directory for output
    afile = args.aln_file #reading in alignment file
    mfile = args.masking_file #reading in masking file
    tree = args.tree #reading in tree name 
    
    print()
    print("Checking dependencies mate! \n")

    if not 'FALSE' in depend_check(dependencies):
        print("I can see all dependencies! \n")
    else:
        print("Mate! Not all required dependencies are loaded.")
        sys.exit()


    ##the project directory that holds the samples
    if args.pdir is not None:
        pdir=args.pdir
    else:
        print()
        print('Input directory path not stated. Stopping Bushfire.')
        print('Need to know their origin story before seeking a destination mate.')
        sys.exit()
    ##the project directory that holds the output
    odir=args.odir

    #reading in reference file
    reference=args.reference


    ##where do you want the output to go

    ##if the project directory and output directory doesn't have a forward slash exit
    if(pdir[-1]!='/'):
      print(pdir[-1])
      print('\n The input directory should end with a forward slash')
      exit()

    ##check if the ids are provided individually or in a file. If they are in a file, read the file and get the ids
    idfile=args.ids[0]
    allids=pd.read_csv(idfile,header=None)
    allids.columns=['ids']
    idlist=allids.ids.tolist()

    #Let your peeps know what is happening. Just a bit of communication.
    print(" ")
    print('Input directory will be: ' + pdir)
    print('Output directory will be: ' + odir)
    print('Using reference file: ' + reference)
    print('Using id csv file: ' + args.ids[0])


    #############################################################################################
    #
    #   Time for the three step program to build a JAR tree
    #
    #############################################################################################

    print("Time to start the Aussie tree cycle of life!")
    
    ######################################################
    # Purge the aln file(tree seed) with fire!  (Run purge module)
    ######################################################
    
    purge()

    ######################################################
    # Time to grow the JAR tree (Run Germie)
    ######################################################        
    
    purgefile = 'masked_core.snp_sites.aln'
    germie()
    
    ######################################################
    # Time to put your tree in a potplant to take away (Run potplant)
    ######################################################  
    
    JARtreefile = str(tree+'.joint.tre')
    potplant()


      
    print("\nAfter the bushfire life will grow again. Did you get a good tree?")

if __name__ == '__main__':
    main()