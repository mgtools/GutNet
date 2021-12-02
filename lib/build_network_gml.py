#!/usr/bin/env python3
"""
builds networkx object from correlation matrix
"""

__author__ = "Tony Lam"
__version__ = "0.1.0"
__license__ = "MIT"


import argparse
import numpy as np
import pandas as pd
import networkx as nx 
import sys, os
import leidenalg 
import igraph as ig

def get_df(infile, ctype, corr_threshold):
    ''' 
    prepares df
    '''
    #filter dataframe based on type and threshold 
    if ctype == 'spiec-easi':
        # load df
        df = pd.read_csv(infile, index_col=0)
        # subset df
        df = df.iloc[:,2:5]
        # reorder df
        df = df[['i_annotated','j_annotated', 'x']]
        # rename columns
        df.columns = ['i_annotated','j_annotated', 'corr']
        # filter
        df = df.loc[abs(df["corr"]) >= corr_threshold,]
    elif ctype == 'flashweave':
        # load df
        df = pd.read_csv(infile, index_col=None, comment='#', sep='\t', header=None)
        # rename columns
        df.columns = ['i_annotated','j_annotated', 'corr']
    else:
        sys.exit('invalid/unsupported type')
    return df.reset_index(drop=True)

def add_node(G, node, annotations):
    '''
    add node with annotations to graph
    '''
    idd = annotations.loc[node,'id']
    kingdom = (annotations.loc[node,'Kingdom']).split('__')[1]
    phylum = (annotations.loc[node,'Phylum']).split('__')[1]
    classs = (annotations.loc[node,'Class']).split('__')[1]
    order = (annotations.loc[node,'Order']).split('__')[1]
    family = (annotations.loc[node,'Family']).split('__')[1]
    genus = (annotations.loc[node,'Genus']).split('__')[1]
    try:
        species = (annotations.loc[node,'Species']).split('__')[1]
    except:
        species = np.nan

    G.add_node(str(idd), name=str(idd), identity=node, Kingdom=kingdom, Phylum=phylum, Class=classs, Order=order, Family=family, Genus=genus, Species=species)

    return G

def add_edge(G, A, B, weight, annotations):
    '''
    add edge + weight
    '''
    idd_A = annotations.loc[A,'id']
    idd_B = annotations.loc[B,'id']
    G.add_edge(str(idd_A),str(idd_B),weight=abs(weight), sign = ('positive' if weight >= 0 else 'negative'))

    return G

def run_leidenalg(G):
    '''
    runs leiden algorithm to find community modules
    '''
    # convert networkx graph to igraph object
    g = ig.Graph.from_networkx(G)

    # leiden alg
    partition = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition)

    return(g, partition)

def add_leiden_attribute(G, partition, g):
    '''
    adds leiden module partition into graph attributes
    '''
    for n, part in enumerate(partition):
        for j in part:
            idd = str(g.vs[j]['name'])
            G.nodes[idd]['leiden_partition'] = n
    return G

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='builds networkx object from correlation output, exports graph file')
    parser.add_argument("-i", "--infile", help='input correlation output file')
    parser.add_argument("-t", "--type", help='specify type of correlation used [options: spiec-easi and flashweave]')
    parser.add_argument("-c", "--correlation", help ='correlation threshold to filter', default=0.1, type=float)
    parser.add_argument("-o", "--output", help = 'output matrix file name')
    #parser.add_argument("-l", "--level", help = 'taxa level abundance table agglomeraged to [options: species,genus]', default='species')
    parser.add_argument("--annotations", help = 'phyloseq annotation dictionary id reference (default: 0_taxa_annotations/phyloseq_id_to_annotation.csv)', default = '0_taxa_annotations/phyloseq_id_to_annotation.csv')
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()

    # load annotation metadata
    annotations = pd.read_csv(args.annotations, index_col=2)

    # get df 
    df = get_df(args.infile, args.type, args.correlation)

    # build network
    G = nx.Graph()

    for index,row in df.iterrows():
        #'id', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        A = row['i_annotated']
        B = row['j_annotated']
        weight = row['corr']
        G = add_node(G, A, annotations)
        G = add_node(G, B, annotations)
        G = add_edge(G, A, B, weight, annotations)

    g, partition = run_leidenalg(G)

    G = add_leiden_attribute(G, partition, g)

    nx.write_gml(G, args.output)
    
