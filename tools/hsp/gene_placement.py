#!/usr/bin/env python2.7
# -*-coding:Utf-8 -*
__author__ = ' Moussa Samb & Maria Bernard  & Geraldine Pascal INRAE - SIGENAE '
__copyright__ = 'Copyright (C) 2020 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs@inrae.fr'
__status__ = 'dev'

#Import
import os
import sys
import argparse
import json
import re
from numpy import median

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# PATH: executable
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH'] 
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) 
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR 
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR
###################
print("Partie1")

#import frogs
from frogsUtils import *
from frogsSequenceIO import * 
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES # hsp.py -i 16S -t placed_seqs.tre -o marker_nsti_predicted.tsv.gz -p 1 -n
                # hsp.py -i EC -t placed_seqs.tre -o EC_predicted.tsv.gz -p 1
                # hsp.py -i KO -t placed_seqs.tre -o KO_predicted.tsv.gz -p 1

#
##################################################################################################################################################
class hsp_16S(Cmd):
    """
    @summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
    @summary: hsp.py program  predict number copie of gene family for each OTU.
    @summary: use 16S, EC and/or KO 
    @see: https://github.com/picrust/picrust2/wiki
    """

    def __init__(self, hsp16Met, in_trait, tree, output, calculate_NSTI, stdout):
       
        Cmd.__init__(self,
                 'hsp.py',
                 'predict gene copy 16S', 
                  "" + hsp_16_Met  +" --in_trait "+ str(in_trait) +" --tree "+ str(tree) +" --output "+ str(output) +" --calculate_NSTI"+str(calculate_NSTI) +' 2> ' + stdout,
                "--version") 
      

    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """

        return Cmd.get_version(self, 'stdout').split()[1].strip()
        

#####################################################   
class hsp_EC(Cmd):
    """
    @summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
    @summary: hsp.py program  predict number copie of gene family for each OTU.
    @summary: use EC 
    @see: https://github.com/picrust/picrust2/wiki
    @commande: hsp.py -i EC -t placed_seqs.tre -o EC_predicted.tsv.gz -p 1
    """
    def __init__(self, hspECMet, in_trait, tree, output, stdout):
       
        Cmd.__init__(self,
                 'hsp.py',
                 'predict gene copy 16S', 
                  "" + hsp_EC_Met  +" --in_trait "+ str(in_trait) +" --tree "+ str(tree) +" --output "+ str(output) +' 2> ' + stdout,
                "--version") 
        

    def get_version(self): 
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """

        return Cmd.get_version(self, 'stdout').split()[1].strip() 

       
#################################################################
class hsp_KO(Cmd):
    """
    @summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
    @summary: hsp.py program  predict number copie of gene family for each OTU.
    @summary: use KO 
    @see: https://github.com/picrust/picrust2/wiki
    @commande: hsp.py -i KO -t placed_seqs.tre -o KO_predicted.tsv.gz -p 1
    """
    def __init__(self, hspKOMet, in_trait, tree, output, stdout):
       
        Cmd.__init__(self,
                 'hsp.py',
                 'predict gene copy 16S', 
                  "" + hsp_KO_Met  +" --in_trait "+ str(in_trait) +" --tree "+ str(tree) +" --output "+ str(output) +' 2> ' + stdout,
                "--version") 
        

    def get_version(self): 
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """

        return Cmd.get_version(self, 'stdout').split()[1].strip() 

  
    ########################
    print("partie2")
##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
    #Fonction restricted float

################################
print("partie3")
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
    # Table to use for prediction
    TRAIT_OPTIONS = ['16S', 'COG', 'EC', 'KO', 'PFAM', 'TIGRFAM', 'PHENO']
    # prediction method to use 
    HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']

    # Manage parameters
    parser = argparse.ArgumentParser( description='predict gene family for OTU' )


    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )

    group_input.add_argument('-i', '--in_trait', type=str.upper, choices=TRAIT_OPTIONS, help='Specifies which default trait table should be used. Use the --observed_trait_table option to input a non-default trait table.')

    group_input.add_argument('-t', '--tree', metavar='PATH', required=True, type=str, help='The full reference tree in newick format containing both study sequences (i.e. ASVs or OTUs) and reference sequences.')

    group_input.add_argument('-o', '--output', metavar='PATH', type=str, required=True, help='Output table with predicted abundances per study sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.')

    group_input.add_argument('--observed_trait_table', metavar='PATH', type=str, help='The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you want to use a custom table.')

    group_input.add_argument('--chunk_size', default=500, type=int, help='Number of functions to run at a time on one processor. Note that you should consider how many ''processes you have specified before changing this ''option. E.g. if you specify the chunk_size to be ''the total number of functions, 1 processor will ''be used even if you specified more so the job will ''be substantially slower (default: %(default)d).')

    group_input.add_argument('-m', '--hsp_method', default='mp',choices=HSP_METHODS, help='HSP method to use.' +'"mp": predict discrete traits using max parsimony. ''"emp_prob": predict discrete traits based on empirical ''state probabilities across tips. "subtree_average": ''predict continuous traits using subtree averaging. ' '"pic": predict continuous traits with phylogentic ' 'independent contrast. "scp": reconstruct continuous ''traits using squared-change parsimony (default: ''%(default)s).')

    group_input.add_argument('-p', '--processes', default=1, type=int, help='Number of processes to run in parallel (default: ' '%(default)d).')

    group_input.add_argument('--verbose', default=False, action='store_true',help='If specified, print out wrapped commands and other ''details to screen.')

    #group_input.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)
  
    # output
    group_output = parser.add_argument_group( 'Outputs' )

    #group_output.add_argument('-m','--html', default='summary.html', help="Path to store resulting html file. [Default: %(default)s]" )    

    group_output.add_argument('-n', '--calculate_NSTI', default=False, action='store_true', help='Calculate NSTI and add to output file.')

    group_output.add_argument('--check', default=False, action='store_true', help='Check input trait table before HSP.')

    group_output.add_argument('--seed', default=100, type=int, help='Seed to make output reproducible, which is ''necessary for the emp_prob method ' '(default: %(default)d).')

    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')


    args = parser.parse_args()
    #######
    stderr = "hsp.stderr"
    ######
    hsp_16_Met = ""
    hsp_EC_Met = ""
    hsp_KO_Met = ""

    # Process 
    try:     
        #~ print(default_ref_dir)
        #~ print(restricted_float)
        print("\n\n\n--------------")
        #~ print(args)
        os.system("pwd")
        #####
        hsp_16S(hsp_16_Met, args.in_trait, args.tree, args.calculate_NSTI, stderr).submit( args.log_file )
        #
        hsp_EC(hsp_EC_Met, args.in_trait, args.tree, stderr).submit( args.log_file )
        #
        hsp_KO(hsp_KO_Met, args.in_trait, args.tree, stderr).submit( args.log_file )
        ########
        print("Partie hsp fini")
        ########
        #write_summary( args.html, "sout.fasta", "excluded.tsv", args.biom_file, args.out_tree)
    finally:
        print("Partie Finale ")

        