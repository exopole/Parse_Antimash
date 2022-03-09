#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
#author          :Annie Lebreton && Alexandre Nicaise
#date            :
#version         :1
#python_version  :3.9
"""

import argparse
import glob
from Bio import SeqIO
from Bio.SeqFeature import (
    FeatureLocation,
)
from collections import Counter


class Gene:
    def __init__(self, name, sequence):
        self.Name = name
        self.Sequence = sequence


class Cluster:
    def __init__(self):
        self.Sp = ""
        self.Tag = ""
        self.location = FeatureLocation(0, 0)
        self.Genes = []

    def ToString(self):
        result = self.GetStringCommun()
        return result

    def GetStringCommun(self):
        result = ""
        result += self.Sp + "\t"
        result += self.Tag + "\t"
        result += str(self.location.start.position) + ":" + str(self.location.end.position) + "\t"
        return result


def parseoptions():
    """ Docstring 
    .... """

    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-d', '--directory')

    global ARGS  # Update the global ARGS variable
    ARGS = parser.parse_args()


##----- INPUT ---------

##----- OUTPUT --------


def get_cds_feature_with_qualifier_value(seq_record, name, value):
    """Function to look for CDS feature by annotation value in sequence record.
    
    e.g. You can use this for finding features by locus tag, gene ID, or protein ID.
    """
    # Loop over the features
    for feature in genome_record.features:
        if feature.type == "CDS" and value in feature.qualifiers.get(name, []):
            return feature
    # Could not find it
    return None


def parse_id_gene(name):
    """
    Funtion to parse a name (jgi.p_Abobie1_7) to get the ID (Abobie1)
    """
    type(name)
    output = name.split("_")
    if len(output) >= 2:
        return output[1]
    return name

def ParsingAndWriting(infile, outfile):
    o = open(outfile, "w")

    clusters = []
    cluster = Cluster()


    with open(infile) as handle:
        for record in SeqIO.parse(handle, format="gb"):
            for feature in record.features:
                if feature.type == "cand_cluster":
                    if cluster.Tag != "" and len(cluster.Genes) > 0:
                        # print("In for => " + str(len(cluster.Genes)))
                        clusters.append(cluster)
                        cluster = Cluster()
                        cluster.Genes = []

                    product = feature.qualifiers["product"][0]
                    cluster.Tag = product
                    cluster.location = feature.location
                if feature.type == "CDS":
                    if "gene_functions" in feature.qualifiers:
                        if product in feature.qualifiers["gene_functions"][0]:
                            if cluster.Sp == "" and "gene" in feature.qualifiers:
                                cluster.Sp = parse_id_gene(feature.qualifiers["gene"][0])
                                # print(feature.qualifiers["gene"])
                                # print(len(cluster.Genes))
                                gene = Gene(str(feature.qualifiers["gene"][0]), str(feature.qualifiers["translation"][0]))
                                cluster.Genes.append(gene)

            if cluster.Tag != "" and len(cluster.Genes) > 0:
                # print("Appen after record ===> " + str(len(cluster.Genes)))
                clusters.append(cluster)
                cluster = Cluster()
                cluster.Genes = []

    toWrite = ""
    for item in clusters:
        for gene in item.Genes:
            toWrite += item.GetStringCommun() + str(gene.Name) + "\n"
            print(">" + str(gene.Name) + "|" + item.Tag + "|" + item.Sp + "|" + str(item.location.start.position) + ":"
                  + str(item.location.end.position))

    o.write(toWrite)
    o.close()


def main():
    parseoptions()
    if ARGS.infile is not None:
        ParsingAndWriting(ARGS.infile, ARGS.outfile)

    if ARGS.directory is not None:
        path = ARGS.directory
        path_to_file = glob.glob(path + '*.gbk')
        for file in path_to_file:
            print(file)
            pathFileSplit = file.split('/')
            nameFile = pathFileSplit[len(pathFileSplit)-1]
            ParsingAndWriting(file, nameFile[:-4])

main()
