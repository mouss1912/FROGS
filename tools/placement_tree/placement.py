#!/usr/bin/python2.7
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
#from cmd import Cmd

################ Ajouter des dossiers dans la variable PATH (libexec),Permet d'éxecuté le script ######################
# repertoire du chemin absolu
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
print(CURRENT_DIR ,"llllll")
# PATH: executable
BIN_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "libexec"))
print(BIN_DIR ,"llllll")
os.environ['PATH'] = BIN_DIR + os.pathsep + os.environ['PATH'] #Permet d'éxécuter (des outils quelques soit l'endroit) #echo $PATH affiche tous les bin
# PYTHONPATH
LIB_DIR = os.path.abspath(os.path.join(os.path.dirname(CURRENT_DIR), "lib"))
sys.path.append(LIB_DIR) #Ajoute ce nouveau chemin dans mon path
if os.getenv('PYTHONPATH') is None: os.environ['PYTHONPATH'] = LIB_DIR #
else: os.environ['PYTHONPATH'] = os.environ['PYTHONPATH'] + os.pathsep + LIB_DIR
#os.environ['PICRUSt2_PATH'] = "" #variable environnement = pour retrouver les variable 
print("Partie1")
#Importer les modules picrust2 (J'ai mis mes librairie picrsut2 dans le dossier lib)
#from picrust2.place_seqs import place_seqs_pipeline
from picrust2.default import default_ref_dir
from picrust2.util import restricted_float


#Ajouter un lien symbolique de place_seq dans libexec (qui )
#which place_seqs.py

from frogsUtils import *
from frogsSequenceIO import * #Pour parser les Sequences(fasta et fastq)
from frogsBiom import BiomIO

##################################################################################################################################################
#
# COMMAND LINES # place_seqs.py --study_fasta --min_align  --out_tree --ref_dir --threads
#
##################################################################################################################################################
class picrust2_place_seqs(Cmd):# la class Picrust2 herite de la class Cmd # Class picrust2_place_seq
    """
    @summary: PICRUSt2 (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States)
    @summary: place_seqs.py program placed Fasta sequence on the tree Arbre .
    @see: https://github.com/picrust/picrust2/wiki
    """
    # methode _init_ = Comme le constructeur de la classe "1er methode à creer", initialise les attributs
    def __init__(self, picrustMet, study_fasta, out_tree,min_align, ref_dir, stdout):
        """
        @param mafftMet: [str] picrust2 method option.
        @param fasta: [str] Path to input fasta file.
        @param out_tree: [str] Path to store resulting tree file.
        @param stderr: [str] Path to temporary picrust2 stderr output file
        @param thread: [int] number of cpu to use.
        """
        # La classe commande n'entre pas dans le programme picrust2
       # print("" + picrustMet  +" --study_fasta "+ str(study_fasta) +" --out_tree "+ str(out_tree) +" --min_align "+ str(min_align) +" --ref_dir "+ str(ref_dir))
        os.system("place_seqs.py "+ picrustMet  +" --study_fasta "+ str(study_fasta) +" --out_tree "+ str(out_tree) +" --min_align "+ str(min_align) +" --ref_dir "+ str(ref_dir))
        """
        Cmd.__init__(self,
                 'place_seqs.py', #l'executable (place_seq.py) 
                 'place OTUs on tree.', # c'est le descriptif de l'outil (place otu on tree)
                  "" + picrustMet  +" --study_fasta "+ str(study_fasta) +" --out_tree "+ str(out_tree) +" --min_align "+ str(min_align) +" --ref_dir "+ str(ref_dir) +' 2> ' + stdout,
                "--version") #Ajouter le ref_dir et le min_align (ajoute version argument)
        """
        """

        print(self.description)
        print(os.path.basename(self.program))
        print(self.get_version())
        """


    
    def get_version(self): #fonction qui fait appel à Cmd et --version
        """
        @summary: Returns the program version number.
        @return: [str] Version number if this is possible, otherwise this method return 'unknown'.
        """

        print (Cmd.get_version(self,'stdout'), "jjkjk")
        print (Cmd.get_version(self,'stdout').split())
        print(Cmd.get_version(self,'stdout')[1])
        return Cmd.get_version(self, 'stdout').split()[1].strip() # Le strip ça enléve le retour chariot # Le stdout prend l'erreur de la sortie standard)

    

##################################################################################################################################################
#
# FUNCTIONS
#
##################################################################################################################################################

def get_fasta_nb_seq( fasta_file ):
    """
    @summary: Returns the number of sequences in fasta_file.
    @param fasta_file: [str] Path to the fasta file to process.
    @return: [int] The number of sequences.
    @ Adapt fasta file to picrust2 input
    """
    return sum(1 for _ in (FastaIO(fasta_file)))



#Changement de la fonction parce que non compatible avec l'entré de Picrust2 (reference et espace pas compris par hmalign)
# Fonction de conversion ddu fichier FASTA
def convert_fasta(fasta_file):

    input_fasta = fasta_file
    output_fasta = "sout.fasta"

    FH_input = FastaIO( input_fasta )
    # pour écrire dans un fichier au format fasta
    f_out = open(output_fasta, "w")
    #FH_output = FastaIO( output_fasta, "w" )
    chaine = ""
    for record in FH_input:
        # récupération du nnom de la séquence 
        #print(record.id)
        # tu peux modifier l'identifiant, par exemple en ajoutant Moussa
        record.id = record.id
        # récupération de la description
        #print(record.description)
        record.description = ""
        f_out = open(output_fasta, "w")
        # récupération de la séquence
        ##print(record.string)
        #print(ide)
        #print(str(ide))
        #print(record.string)
        chaine += ">"+record.id+"\n"+record.string+"\n"
        # pour écrire la sequence (identifiant, description, et sequence) dans le fichier de sortie)
        
        #FH_output.write(record)
    f_out.write(chaine)
    f_out.close()
    



#Fonction pour Extraire les cluster non alignés

def excluded_sequence(file_tree, file_fasta, out_file):

    """
    @summary: Returns the excluded sequence.
    @param fasta_file: [str] Path to the fasta file to process.
    @param tree_file: [str] Path to the tree file to process.
    @return: [int] The file of no aligned sequence.
    """
    #Lecture du fichier tree
    file = open(file_tree, "r")
    line = file.readline()
    #List of cluster
    list_cluster = []
    #Boucle sur le fichier tree
    #Je splite sur (,) , regex sur le cluster et récupérer le groupe1
    while line: 
        for i, v in enumerate(line.split(",")):
            group = re.search("(Cluster_[0-9]+)", v)
            if group:
                ide = group.group(1)
                list_cluster.append(ide)
        line = file.readline() #parcourir ligne par ligne
    file.close()

    file_fasta = open(file_fasta, "r")
    file_out   = open(out_file, "w")

    line = file_fasta.readline()
    while line:
        if line[0] == ">":
            if line[1:].strip() not in list_cluster:
                file_out.write(line[1:].strip()+"\n")

        line = file_fasta.readline()

    file_fasta.close()
    file_out.close()

#Fonction pour le summary html
#def write_summary( summary_file, file_fasta, out_tree, file_tree):
    """
    @summary: Writes the process summary in one html file.
    @param summary_file: [str] path to the output html file.
    @param file_fasta: [str] path to the fasta file of OTU
    @param out_tree [str] path to the input tree file.
    @param file_tree: [str] path to the Newick file.
    """
##################################################################################################################################################
#
# MAIN
#
##################################################################################################################################################

if __name__ == "__main__":

    # Manage parameters
    parser = argparse.ArgumentParser( description='Phylogenetic tree reconstruction' )
    # Inputs
    group_input = parser.add_argument_group( 'Inputs' )

    group_input.add_argument('-s', '--study_fasta', metavar='PATH', required=True,type=str, help='FASTA of unaligned study sequences.')

    group_input.add_argument('-p', '--processes', type=int, default=1, help='Number of processes to run in parallel (default: ''%(default)d). Note that this refers to ''multithreading rather than multiprocessing when ''running EPA-ng and GAPPA.')

    group_input.add_argument('-r', '--ref_dir', metavar='PATH', type=str, default=default_ref_dir, help='Directory containing reference sequence files ''(default: %(default)s). Please see the online ''documentation for how to name the files in this ''directory in order to use custom reference files.')

    group_input.add_argument('--min_align', type=restricted_float, default=0.8, help='Proportion of the total length of an input query ''sequence that must align with reference sequences. ''Any sequences with lengths below this value after ''making an alignment with reference sequences will ''be excluded from the placement and all subsequent ''steps. (default: %(default)d).')

    group_input.add_argument('--chunk_size', type=int, default=5000, help='Number of query seqs to read in at once for EPA-ng ''(default: %(default)d).')

    group_input.add_argument('--verbose', default=False, action='store_true', help='If specified, print out wrapped commands and other ''details to screen.')

    group_input.add_argument('-v', '--version', default=False, action='version', version="%(prog)s " + __version__)

  
    # output
    group_output = parser.add_argument_group( 'Outputs' )

    group_output.add_argument('-o', '--out_tree', metavar='PATH', required=True, type=str, help='Name of final output tree.')

    group_output.add_argument('--intermediate', metavar='PATH', type=str, default=None, help='Output folder for intermediate files (will be ''deleted otherwise).')
    #group_output.add_argument( '-l', '--log-file', default=sys.stdout, help='This output file will contain several information on executed commands.')

    args = parser.parse_args()
    #prevent_shell_injections(args)

     ### Temporary files
    # alignment temporary files
    #Pas bon pas " code programme doit fonctiooné partout"
    stderr = "picrust2.stderr"

    picrustMet = ""
    # Process 
    try:     
        print(default_ref_dir)
        print(restricted_float)
        print("\n\n\n--------------")
        print(args)
        os.system("pwd")
        #place_seqs.py --study_fasta --min_align  --out_tree --ref_dir --threads
        convert_fasta(args.study_fasta)
        picrust2_place_seqs(picrustMet, "sout.fasta", args.out_tree, args.min_align, args.ref_dir, stderr)
        print("Partie picrust fini")
        excluded_sequence(args.out_tree, "sout.fasta", "excluded.tsv")
        #PICRUSt2(picrustMet, "sout.fasta", args.out_tree, args.min_align, args.ref_dir, stderr).submit( args.log_file )
        print("Partie 2 ")
    finally:
    	print("Partie Finale ")
        

    