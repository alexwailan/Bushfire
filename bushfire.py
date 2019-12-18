#!/usr/bin/env python3


#Created: 27.11.18 - Alexander Wailan 



import os
import sys
import gzip
import argparse
import pandas as pd
import numpy as np
import csv
import subprocess
import pathlib
from pathlib import Path
import shutil

##########################################
# Check if dependent programs are loaded #
##########################################

dependencies = [
'samtools',
'bcftools'
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

print()
print("Checking dependencies mate! \n")

if not 'FALSE' in depend_check(dependencies):
    print("I can see all dependencies! \n")
else:
    print("Mate! Not all required dependencies are loaded.")
    sys.exit()


##########################################
# Function to Get command line arguments #
##########################################

def getargv():
    usage = 'bushbush.py [options] reference ids'
    description='Run Bushfire. A program to parse a multi-fasta core SNV alignment file to a tree'
    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('reference', action="store", help='Provide the reference for snippy [Required]', metavar='N',nargs='?')
    parser.add_argument('ids', action="store", help='IDs of paired reads files [Required]', metavar='N', type=str, nargs='+')
    parser.add_argument('-d',    '--dirpath',  action="store",dest="pdir", help='Input directory of Lodestone results. End with a forward slash. Eg. /temp/fasta/ [Required]', metavar='N', type=str, nargs='?')
    parser.add_argument('-o',    '--outdir', action="store", dest="odir", help='Output directory. End with a forward slash. Eg. /temp/fasta/; Default to use current directory.', metavar='N', type=str, nargs='?', default=os.getcwd()+'/')
    return parser.parse_args()


args = getargv()

################
# Main program #
################


#############################################################################################
#
#      Parse/ check the arguements
#
#############################################################################################


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
#      Construct the output command for the program, use the slurm script as a template
#      to construct the shell script & run
#
#############################################################################################

print("Time to start the Aussie tree cycle of life!")
for i in idlist:
    print(" ")
    print("Bushfire with id: %s"%i)
    if not (i in os.listdir(odir)):
        os.system(""" mkdir %s"""%i) #directory structure required for snippy-core
        print("Creating directory for isolate %s" %i)

    if (i in os.listdir(odir)):
        
        ######################################################
        # copy lodestone vcf output into each sample directory
        ######################################################
        
        #WILL HAVE TO CHANGE THIS DIRECTORY
        
        p = subprocess.call("cp -r %s %s"%(pdir+i+'.vcf.gz', odir+i+'/'), shell=True,    stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        gzip_f = Path(odir+i+'/'+i+'.vcf.gz')

        if not os.path.isfile(gzip_f):
            print("Stopped Bushwalking with %s." %i)
            print("Copying files failed. Unable to find required vcf.gz file for bcftools." %i)
            break
        
        ######################################################
        # Create VCF file with only SNVs
        ######################################################        

        p = subprocess.Popen("bcftools view --types snps %s > %s"%(odir+i+'/'+i+'.vcf.gz',odir+i+'/'+i+'.SNV.vcf'), shell=True,    stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        output,error = p.communicate() #Read data from stdout and stderr. Wait for process to terminate.
        vcf_f = Path(odir+i+'/'+i+'.SNV.vcf')
        
        if not os.path.isfile(vcf_f):
            print("Stopped Bushwalking with %s." %i)
            print("bcftools view failed")
            print(error)
            break
        
        ######################################################
        # Compress for indexing & Index vcf
        ######################################################   

        p = subprocess.Popen("bgzip %s"%(odir+i+'/'+i+'.SNV.vcf'), shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        output,error = p.communicate() #Read data from stdout and stderr. Wait for process to terminate.
        bg_f = Path(odir+i+'/'+i+'.SNV.vcf.gz')
        
        if not os.path.isfile(bg_f):
            print("Stopped Bushwalking with %s." %i)
            print("bgzip failed")
            break

            #index the vcf file
        p = subprocess.Popen("tabix -p vcf %s"%(odir+i+'/'+i+'.SNV.vcf.gz'), shell=True,    stdout=subprocess.PIPE,stderr=subprocess.PIPE)

        
        index_f = Path(odir+i+'/'+i+'.SNV.vcf.gz.tbi')
        output,error = p.communicate() #Read data from stdout and stderr. Wait for process to terminate.
        
        if not os.path.isfile(index_f):
            print("Stopped Bushwalking with %s." %i)
            print("tabix failed")
            print(error)
            break

        ######################################################
        # Generate consensus alignment file with bcftools
        ######################################################   

        p = subprocess.Popen("bcftools consensus -f %s -o %s %s"%(reference,odir+i+'/'+i+'.SNV.aligned.fa',odir+i+'/'+i+'.SNV.vcf.gz'), shell=True,    stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        #print(odir+i+'/'+i+'.SNV.vcf.gz')
        output,error = p.communicate() #Read data from stdout and stderr. Wait for process to terminate.
        aln_f = Path(odir+i+'/'+i+'.SNV.aligned.fa')

        if not os.path.isfile(aln_f):
            print("Stopped Bushwalking with %s." %i)
            print("concensus failed")
            print(error)
            break

        ######################################################
        # Decompress vcf.gz to be used for snippy-core
        ######################################################   

        p = subprocess.Popen("bgzip -d %s "%(odir+i+'/'+i+'.SNV.vcf.gz'), shell=True,    stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        output,error = p.communicate() #Read data from stdout and stderr. Wait for process to terminate.
        vcf_f = Path(odir+i+'/'+i+'.SNV.vcf')

        if not os.path.isfile(vcf_f):
            print("Stopped Bushwalking with %s." %i)
            print("bgzip decompression failed")
            print(error)
            break
        p = subprocess.Popen("rm  %s "%(odir+i+'/'+i+'.vcf.gz'), shell=True,    stdout=subprocess.PIPE,stderr=subprocess.PIPE) #Clean up
        
        if os.path.isfile(vcf_f) and os.path.isfile(aln_f):
            print("Bushwalk was successfull for %s. "%i)
print()
print("Bushwalk complete. Did you find the trees?")
