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

#import argparse
#from picrust2.wrap_hsp import castor_hsp_workflow
#from picrust2.util import make_output_dir_for_file, check_files_exist
#from picrust2.default import default_tables

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
#os.system("hsp.py "+ hsp_16_Met  +" --in_trait "+ str(in_trait) +" --tree "+ str(tree) +" --output "+ str(output) +" --calculate_NSTI"+str(calculate_NSTI))

class hsp16S(Cmd):
    """
    @summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
    @summary: hsp.py program  predict number copie of gene family for each OTU.
    @summary: use 16S, EC and/or KO 
    @see: https://github.com/picrust/picrust2/wiki
    """
    
    def __init__(self, hsp16Met,categorie, in_trait, tree, output, n, stdout):
    #def __init__(self, hsp16Met, categorie, in_trait, tree, output, n, stdout):    
       #os.system("hsp.py "+ hsp16Met  +" -i "+ str(in_trait) +" -t "+ str(tree) +" -o "+ str(output) +" --observed_trait_table "+ str(observed_trait_table) + " -n " )
        #-n: Calculate NSTI and add to output file.
        Cmd.__init__(self,
                 'hsp.py',
                 'predict gene copy 16S', 
                  "" + hsp16SMet  +" -i "+ str(in_trait) +" -t "+ str(tree) +" -o "+ str(output) +" -n" +' 2> ' + stdout,
                "--version") 
      
        
    def get_version(self):
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """

        return Cmd.get_version(self, 'stdout').split()[1].strip()
        

#####################################################   
class hspITS(Cmd):
    """
    @summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
    @summary: hsp.py program  predict number copie of gene family for each OTU.
    @summary: use EC 
    @see: https://github.com/picrust/picrust2/wiki
    @commande: hsp.py -i EC -t placed_seqs.tre -o EC_predicted.tsv.gz -p 1
    """
    def __init__(self, hspITSMet,categorie, tree, output, observed_trait_table, n, stdout):
        #os.system("hsp.py "+ hsp_EC_Met  +" -i "+ str(in_trait) +" -t "+ str(tree) +" -o "+ str(output))
        #os.system("hsp.py "+ hspITSMet+ " -t "+ str(tree) +" -o "+ str(output) +" --observed_trait_table "+ str(observed) + " -n ")

       
        Cmd.__init__(self,
                 'hsp.py',
                 'predict gene copy fungi', 
                  "" + hspITSMet  +" -t "+ str(tree) +" -o "+ str(output) + " --observed_trait_table "+ str(observed_trait_table) +" -n" +' 2> ' + stdout,
                "--version") 
    
    def get_version(self): 
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """

        return Cmd.get_version(self, 'stdout').split()[1].strip() 



##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################
def write_summary( summary_file, output, biom_file):
    """
    @summary: Writes the summary of results.
    @param summary_file: [str] The output file.
    @param results_chimera: [str] Path to the input chimera step summary.
    """
    # Get data
    detection_categories = ["sequence", "16S_rRNA_Count"]
    detection_data = list()
    remove_data = dict()
    """
    # Parse results chimera
    in_remove_metrics = True
    in_detection_metrics = False
    section_first_line = True
    log_fh = open(output)
    
    for line in log_fh:
        line = line.strip()
        if line.startswith('##Metrics by sample'):
            remove_metrics = False
            in_detection_metrics = True
            section_first_line = True
        elif line.startswith('##Metrics global'):
            remove_metrics = True
            in_detection_metrics = False
            section_first_line = True
        elif line == "":
            in_detection_metrics = False
            in_remove_metrics = False
        else:
            if in_detection_metrics:
                if section_first_line:
                    line_fields = line[1:].split("\t")[1:]
                    detection_categories = line_fields
                    section_first_line = False
                else:
                    line_fields = line.split("\t")
                    detection_data.append({
                             'name': line_fields[0],
                             'data': map(int, line_fields[1:])
                    })
            elif in_remove_metrics:
                if section_first_line:
                    line_fields = line[1:].split("\t")
                    remove_categories = [category.lower().replace(" ", "_") for category in line_fields]
                    section_first_line = False
                else:
                    for idx, val in enumerate(line.split("\t")):
                        remove_data[remove_categories[idx]] = int(val)
    
    log_fh.close()
    """
    # Write
    FH_summary_tpl = open( os.path.join(CURRENT_DIR, "gene_placement.html") )
    FH_summary_out = open( summary_file, "w" )
    for line in FH_summary_tpl:
        if "###DETECTION_CATEGORIES###" in line:
            line = line.replace( "###DETECTION_CATEGORIES###", json.dumps(detection_categories) )
        elif "###DETECTION_DATA###" in line:
            line = line.replace( "###DETECTION_DATA###", json.dumps(detection_data) )
        elif "###REMOVE_DATA###" in line:
            line = line.replace( "###REMOVE_DATA###", json.dumps(remove_data) )
        print(line)
        FH_summary_out.write( line )

    FH_summary_out.close()
    FH_summary_tpl.close()
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":
    #Categorie
    Categorie = ['16S', 'ITS', '18S']
    # Table to use for prediction
    TRAIT_OPTIONS = ['16S', 'COG', 'EC', 'KO', 'PFAM', 'TIGRFAM', 'PHENO']
    # prediction method to use 
    HSP_METHODS = ['mp', 'emp_prob', 'pic', 'scp', 'subtree_average']
    #Pour personaliser m
    #retourne le couple de valeur: 1er indice , 2em la valeur
    find = False
    for i, v in enumerate(sys.argv):
        if sys.argv[i] == "-c" or sys.argv[i] == "--categorie":
            if sys.argv[i+1] == "16S":
                find = True      

    # Manage parameters
    parser = argparse.ArgumentParser( description='predict gene family for OTU' )

    print("partie4")
    # Inputs
    #Les inputs categoriel
    group_input = parser.add_argument_group( 'Inputs' )

    group_input.add_argument('-c', '--categorie',choices=Categorie, help='Specifies which categorie 16S or ITS, 18S')

    if find:
        group_input.add_argument('-i', '--in_trait',choices=TRAIT_OPTIONS,  help='Specifies which default trait table should be used. Use the --observed_trait_table option to input a non-default trait table.')

    else:
        group_input.add_argument('-i', '--observed_trait_table', metavar='PATH', type=str, help='The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you want to use a custom table.')

    #Pour avoir plusieur parametres   
    # tab_arg = []
    # for i, v in enumerate(sys.argv):
    #      if sys.argv[i] == "-i":
    #         ind = 0
    #         for e in range(i+1, len(sys.argv)):
    #             if sys.argv[e][0] == "-" :
    #                 ind = e
    #                 print(i+1, ind, len(sys.argv), sys.argv[e])
    #                 break

    #         if ind != 0:
    #             tab_arg = sys.argv[i+1:ind]
    #             print(tab_arg)
    #             break

    #group_input.add_argument('-t', '--table', help='Specifies which table choices, ')

    #les parametres intermediares
    #group_input.add_argument('-i', '--in_trait', type=str.upper, choices=TRAIT_OPTIONS, help='Specifies which default trait table should be used. Use the --observed_trait_table option to input a non-default trait table.')
    #group_input.add_argument('-i', '--in_trait',  help='Specifies which default trait table should be used. Use the --observed_trait_table option to input a non-default trait table.')

    group_input.add_argument('-t', '--tree', metavar='PATH', required=True, type=str, help='The full reference tree in newick format containing both study sequences (i.e. ASVs or OTUs) and reference sequences.')

    group_input.add_argument('-o', '--output', metavar='PATH', type=str, required=True, help='Output table with predicted abundances per study sequence in input tree. If the extension \".gz\" is added the table will automatically be gzipped.')

    #group_input.add_argument('--observed', metavar='PATH', type=str, help='The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you want to use a custom table.')
    #group_input.add_argument('--observed_trait_table', metavar='PATH', type=str, help='The input trait table describing directly observed traits (e.g. sequenced genomes) in tab-delimited format. Necessary if you want to use a custom table.')


    group_input.add_argument('--chunk_size', default=500, type=int, help='Number of functions to run at a time on one processor. Note that you should consider how many ''processes you have specified before changing this ''option. E.g. if you specify the chunk_size to be ''the total number of functions, 1 processor will ''be used even if you specified more so the job will ''be substantially slower (default: %(default)d).')

    group_input.add_argument('-m', '--hsp_method', default='mp',choices=HSP_METHODS, help='HSP method to use.' +'"mp": predict discrete traits using max parsimony. ''"emp_prob": predict discrete traits based on empirical ''state probabilities across tips. "subtree_average": ''predict continuous traits using subtree averaging. ' '"pic": predict continuous traits with phylogentic ' 'independent contrast. "scp": reconstruct continuous ''traits using squared-change parsimony (default: ''%(default)s).')

    group_input.add_argument('-p', '--processes', default=1, type=int, help='Number of processes to run in parallel (default: ' '%(default)d).')

    group_input.add_argument('--verbose', default=False, action='store_true',help='If specified, print out wrapped commands and other ''details to screen.')

    group_input.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)

    group_input.add_argument('-b', '--biom_file', metavar='PATH', required=True, type=str, help='Biom file.')

  
    # output
    group_output = parser.add_argument_group( 'Outputs' )

    group_output.add_argument('-s','--html', default='summary.html', help="Path to store resulting html file. [Default: %(default)s]" )    

    group_output.add_argument('-n', '--calculate', default=False, action='store_true', help='Calculate NSTI and add to output file.')

    group_output.add_argument('--check', default=False, action='store_true', help='Check input trait table before HSP.')

    group_output.add_argument('--seed', default=100, type=int, help='Seed to make output reproducible, which is ''necessary for the emp_prob method ' '(default: %(default)d).')

    group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')

    print("partie5")
    args = parser.parse_args()
    tab_arg = []
    if find:
        tab_arg = args.in_trait.split(",")
    else :
        tab_arg = args.observed_trait_table.split(",")

    print(tab_arg)
    ####Declarer la sortie d'erreur ###
    stderr = "hsp.stderr"
    ### Declarer les methodes hsp ###
    hsp16SMet = ""
    hspITSMet = ""
    print("partie6")
    # Process 
    try:     
        Logger.static_write(args.log_file, "## Application\nSoftware :" + sys.argv[0] + " (version : " + str(__version__) + ")\nCommand : " + " ".join(sys.argv) + "\n\n")
        # Commands execution
        print("partie7")
        for i, v in enumerate(tab_arg):
            if args.categorie == "16S":
                    hsp16_cmd = hsp16S(hsp16SMet ,"16S", v, args.tree, v+"_"+args.output, args.calculate, stderr).submit(args.log_file)
            else :
                    hspITS_cmd = hspITS(hspITSMet, args.categorie, args.tree, "fichier_"+str(i)+"_"+args.output, v, args.calculate, stderr).submit(args.log_file)         
              
    
        #hsp_16S(hsp_16_Met, args.in_trait, args.tree, args.output, args.observed_trait_table, args.calculate, stderr)
        #hsp16S(hsp16SMet, args.in_trait, args.tree,args.output,args.calculate, stderr).submit(args.log_file)
        #
        print("partie8")
        #hsp_EC(hsp_EC_Met, args.tree, args.output, args.observed_trait_table, args.calculate, stderr)

        
        #hspITS(hspITSMet, args.tree, args.output,  args.observed_trait_table, args.calculate, stderr).submit(args.log_file)
        #hspITS(hspITSMet, args.tree, args.output,  args.observed, args.calculate, stderr)

        print("partie9")    
        write_summary(args.html, "output", args.biom_file)

    finally:
        print("Partie Finale ")

        