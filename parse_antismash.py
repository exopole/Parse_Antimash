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
import sys

from Bio import SeqIO
from Bio.SeqFeature import (
    FeatureLocation,
)
from collections import Counter

geneFunctions = []
geneQualifier = []


def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()  # As suggested by Rom Ruben


class Gene:
    Location: FeatureLocation

    def __init__(self, name, sequence, function, location):
        self.Name = name
        self.Sequence = sequence
        self.Function = function
        self.Location = location

    def tostring(self):
        return self.Name + " " + str(self.Function) + " " + str(self.Location.start.position) + " : " + str(
            self.Location.end.position)


class Cluster:
    def __init__(self):
        self.Sp = ""
        self.Tag = ""
        self.Scaffold = ""
        self.location = FeatureLocation(0, 0)
        self.Genes = []

    def ToString(self):
        result = self.GetStringCommun()
        return result

    def GetStringCommun(self):
        result = ""
        result += self.Sp + "\t"
        result += self.Tag + "\t"
        result += self.Scaffold + "\t"
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


def find_clusters_for_gene(gene, clusters):
    """

    @type clusters: Cluster[]
    @type gene: Gene
    """
    i = 0
    cluster = None
    while i < len(clusters):
        cluster = clusters[i]
        # print(cluster.location.start + " <= " + gene.Location.start + " and " + cluster.location.end + " >= " + gene.Location.end)
        if (cluster.location.start <= gene.Location.start and cluster.location.end >= gene.Location.end):
            cluster.Genes.append(gene)
            if cluster.Sp == "":
                cluster.Sp = parse_id_gene(gene.Name)
        i += 1


def ParsingAndWriting(infile, outfile):
    o = open(outfile, "w")
    # print("files ==> " + infile + " / " + outfile)

    clusters = []
    genes = []
    scaffold = ""

    with open(infile) as handle:
        for record in SeqIO.parse(handle, format="gb"):
            scaffold = record.name
            #print(scaffold)

            for feature in record.features:
                if feature.type == "cand_cluster":
                    # print("In for => " + str(len(cluster.Genes)))
                    cluster = Cluster()
                    cluster.Scaffold = scaffold
                    cluster.Genes = []

                    # print("Location cluster : " + str(feature.location.start))

                    product = feature.qualifiers["product"][0]
                    # print(product)
                    cluster.Tag = product
                    cluster.location = feature.location
                    clusters.append(cluster)

                if feature.type == "CDS":
                    for qualifier in feature.qualifiers:
                        if qualifier not in geneQualifier:
                            geneQualifier.append(qualifier)

                    name = str(feature.qualifiers["gene"][0])
                    sequence = str(feature.qualifiers["translation"][0])
                    function = ""

                    if "gene_functions" in feature.qualifiers:
                        function = feature.qualifiers["gene_functions"]
                    elif "NRPS_PKS" in feature.qualifiers:
                        function = feature.qualifiers["NRPS_PKS"][len(feature.qualifiers["NRPS_PKS"]) - 1].replace(
                            "type:", "").replace(" ", "")
                        # print(function)

                    gene = Gene(name, sequence, function, feature.location)
                    genes.append(gene)

    for gene in genes:
        find_clusters_for_gene(gene, clusters)

    toWrite = ""
    clusterCount = 0
    geneCount = 0
    for item in clusters:
        # print(str(len(item.Genes)))
        for gene in item.Genes:
            for function in gene.Function:
                if item.Tag in function:
                    toWrite += item.GetStringCommun() + str(gene.Name) + "\n"
                    break
            geneCount += 1
            # print(">" + str(gene.Name) + "|" + item.Tag + "|" + item.Sp + "|" + str(item.location.start.position) + ":"
            #       + str(item.location.end.position))
        clusterCount += 1

    o.write(toWrite)
    o.close()


def main():
    parseoptions()
    if ARGS.infile is not None:
        ParsingAndWriting(ARGS.infile, ARGS.outfile)

    if ARGS.directory is not None:
        path = ARGS.directory
        path_to_file = glob.glob(path + '/*.gbk')
        count = 0
        progress(count, len(path_to_file), path_to_file[0])

        for file in path_to_file:
            pathFileSplit = file.split('/')
            nameFile = pathFileSplit[len(pathFileSplit) - 1]
            ParsingAndWriting(file, ARGS.outfile + nameFile[:-4])
            count += 1
            progress(count, len(path_to_file), file)

    # for qualifier in geneQualifier:
    #     print(qualifier)

    # for function in geneFunctions:
    #     print(function)


main()
