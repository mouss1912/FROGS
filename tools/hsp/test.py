import os
import sys
import argparse
import json
import re

#Dico chaque nom de colonne est une clé: Ma premiére ligne
#Je parcours le fichier je 
dico = {}

def result_file(output):

    file  = open(output, "r")
    line = file.readline()

    global dico
    tab_max = []
    maxe = 0
    if len(dico):
        for i in dico:
            tab_max.append(len(dico[i]))

    if len(tab_max):
        maxe = max(tab_max)

    print("entete" +"_"+ line)
    tab_key = line.strip().split("\t")
    tab_col = []
    dico_2 = {}
    for i in tab_key:
        if i not in dico and maxe != 0:
            dico[i] = [":"] * maxe
        elif i not in dico:
            dico[i] = []

            print("Voila le tableau")
        dico_2[i] = i
        #print(dico)

    for i in dico:
        if i in dico_2:
            tab_col.append(i)


    print("tab col : ", len(tab_col))
    for line in file:
        tab = line.strip().split("\t")
        print(line)
        e = 0
        dep  = 0
        row  = 0
        done = False
        while not done:

            for key in dico:
                print("e :", e)
                if e < len(tab_col) and  key == tab_col[e]:
                    if e < len(tab):
                        dico[tab_col[e]].append(tab[e])
                    else:
                        dico[tab_col[e]].append(":")
                    e += 1
                else:
                    dico[key].append(":")

            if row >= len(dico[key]):
                done = True

            if dep < maxe:
                dep += 1

            row += 1

    file.close()


for output in ["test1.txt", "test2.txt"]:
    result_file(output)


fout = open("Mousse.txt", "w")
fout.write( "\t".join(list(dico.keys())) + "\n" )

row = 0
p = list(dico.keys())[0]
size = len(dico[p])

while row < size :
    col = 0
    for k in dico:
        if col < len(dico[k]):
            fout.write(dico[k][col] + "\t")
        else:
            fout.write(":" + "\t")
        col += 1

    fout.write("\n")
    row += 1

fout.close()

# ft = open("dicoo.txt", "w")
# ft.write(str(dico))
# ft.close()

