#!/usr/bin/python

##Spencer Krieger

import os
import os.path
import sys
from optparse import OptionParser
import pickle as pkl

from pathheuristic import *
from enumeratepaths import *

from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms.directed_paths import b_visit


## GLOBAL VARIABLES
ROOTDIR = ''


def main(args):
    '''
    Makes the hypergraph object, runs the hyperpath heuristic,
    '''
    opts = parseOptions(args)

    H,_,__ = make_hypergraph(ROOTDIR+'/parsed/{}'.format(opts.name))
    H_sources,H_targets,H_high_penalty_sources,H_name_dict = getSourcesTargets(opts.name,H,'hypergraph',opts.source,opts.target)
    source,target = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name)

    # run
    ilpname = '%s/heuristic.lp' % (ROOTDIR,opts.name,opts.type)
    outprefix = '%s/heuristic' % (ROOTDIR,opts.name,opts.type)

    if opts.heuristic:

        taildistancelist,H2,P = tail_path_heuristic(H,source,target)
        pkl.dump(taildistancelist,open('taildistancelist.pkl','wb'))

    if opts.enumeration:
        taildistancelist,H2,P = tail_path_heuristic(H,source,target)
        p,bestp,bestweight = enumerateallpaths(H2,P,source,target)
        print('enumeration results')
        print(p)
        print(bestp)
        print(bestweight)


#############################################

def parseOptions(args):
    '''
    Parses the arguments.
    Returns the parsed arguments list
    '''
    desc = 'python master-script.py [options]'
    parser = OptionParser(usage=desc)

    # General Options
    parser.add_option('','--force',action='store_true',help='Overwrite files if they exist.')
    parser.add_option('','--printonly',action='store_true',help='Print commands to screen, but do not execute them.')
    
    # EXPERIMENTS/TESTS
    parser.add_option('','--enumeration',action='store_true',help='Enumerate all possible source-target paths')
    parser.add_option('','--heuristic',action='store_true',help='Compute heuristic shortest possibly-cyclic hyperpath from S to T.')
    parser.add_option('','--name',type='string',default='WNT5A',help='Name of dataset (WNT5A, CTNNB1, WNT, or ALL). Default=WNT5A.')
    parser.add_option('','--source',type='string',action='append',help='Sources. Default = WNT5A')
    parser.add_option('','--target',type='string',action='append',help='Targets. Default = CTNNB1')

    opts,args = parser.parse_args()

    if not opts.source:
       opts.source = ['WNT5A']
    if not opts.target:
       opts.target = ['CTNNB1']
    print 'OPTIONS ARE',opts
    return opts


def make_hypergraph(file_prefix,delim=';',sep='\t',keep_singleton_nodes=False):
    '''
    Creates the directed hypergraph object from the hypernodes and hyperedges files.
    Returns the hypergraph object.
    '''
    #NOTE: This is directly copied from https://github.com/annaritz/pathway-connectivity
    hypernodes = {}
    with open(file_prefix+'-hypernodes.txt') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split(sep)
            if len(row) == 1:
                hypernodes[row[0]] = ['OtherComplexes-FIX']
            else:
                hypernodes[row[0]] = row[1].split(delim)
    print('%d hypernodes from hypernodes file' % (len(hypernodes)))
    identifier2id = {}
    id2identifier = {}
    H = DirectedHypergraph()
    if keep_singleton_nodes:
        for n in hypernodes:
            H.add_node(n)

    skipped1 = 0
    skipped2 = 0
    tailsizes = []
    headsizes = []
    selfloops = []
    noinselfloops = 0
    indegree = []
    outdegree = []
    numtargets = 0
    with open(file_prefix+'-hyperedges.txt') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split(sep)
            tail = set()
            head = set()

            ## Tail includes tail and regulators.
            ## Head includes head.
            if row[0] != 'None':
                tail.update(row[0].split(delim))
            if row[1] != 'None':
                head.update(row[1].split(delim))
            if row[2] != 'None':
                tail.update(row[2].split(delim))
            #NOTE: This line was commented out to exclude negative regulators in the tail
            #if row[3] != 'None':
                #tail.update(row[3].split(delim))
            hedge_id = row[4]

            ## THIS IS A HACK FOR NOW ( should be incorporated in the make-hypergraph.py code)
            ## IGnore any reactions that have a Reactome Identifier (e.g. has "HSA") instead of
            ## a PAthway Commons identifier.
            if any(['HSA' in s for s in tail]+['HSA' in s for s in head]):
                skipped1+=1
            elif len(tail)==0 or len(head)==0:
                skipped2+=1
            else:
                hid = H.add_hyperedge(tail,head,identifier=hedge_id)
                tailsizes.append(len(tail))
                headsizes.append(len(head))
                intersection = tail.intersection(head)
                if len(intersection) > 0:
                    selfloops.append([v for v in intersection])

                identifier2id[hedge_id] = hid
                id2identifier[hid] = hedge_id

    print('%d reactions skipped because of Reactome identifier' % (skipped1))
    print('%d reactions skipped because of an empty tail or head' % (skipped2))
    ## annotate nodes
    num_hypernodes = 0
    for node in H.get_node_set():
        if node in hypernodes and hypernodes[node] != [node]:
            H.add_node(node,hypernode_members=hypernodes[node],is_hypernode=True)
            num_hypernodes+=1
        else:
            H.add_node(node,is_hypernode=False,hypernode_members=[])

        H.add_node(node)

    return H, identifier2id, id2identifier


def getSourcesTargets(name,H,graph_type,source_list,target_list):
    '''
    sets the sources and targets that were given on the command line
    returns the sources and targets
    '''
    if name == 'WNT5A':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee'])
        #sources = set(['http://pathwaycommons.org/pc12/Protein_037001d14ad7601b82c325eeaac1cc36']) #example 2 to break the loop
        #sources = set(['http://pathwaycommons.org/pc12/Protein_bdf326eedb65f7fe6a0259c5cb8c4ed4','http://pathwaycommons.org/pc12/Protein_037001d14ad7601b82c325eeaac1cc36']) #Trying to find cyclic

        targets = set(['http://pathwaycommons.org/pc12/Protein_6eae9e1fb8e906b20a3ebdf4485a4a3d']) #example 1
        #targets = set(['http://pathwaycommons.org/pc12/Protein_bc47c96d7c652d22f94260b30d5c8043']) #example 2
        #targets = set(['http://pathwaycommons.org/pc12/Complex_46f99a13ca1b39c8d93926b9b394c395'])# trying to find cyclic
    elif name == 'WNT':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee','http://pathwaycommons.org/pc12/Protein_355a3029f445775a6c82451d5c86031b','http://pathwaycommons.org/pc12/Protein_6b903a8964fcd5874a866c1e38e37456'])
        targets = set(['http://pathwaycommons.org/pc12/Protein_3c9a7ce5eec5c6ec97c1a3008c3c9c99','http://pathwaycommons.org/pc12/Protein_1625818ba862a465e6bfe45c1a57c0ec'])
    elif name == 'allpid':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee','http://pathwaycommons.org/pc12/Protein_355a3029f445775a6c82451d5c86031b','http://pathwaycommons.org/pc12/Protein_6b903a8964fcd5874a866c1e38e37456'])
        targets = set(['http://pathwaycommons.org/pc12/Complex_81ba3b0707b6c6abd477dd59315147f4'])
        if ':' in target_list[0]:
            targets = set([target_list[0]])

    print ('sources,targets',sources,targets)

    ## backward star sources
    high_penalty_sources = set([n for n in H.get_node_set() if len(H.get_backward_star(n))==0]).difference(sources)
    all_sources = sources 
    if len(all_sources.intersection(targets)) >0:
        print 'Warning: removing %d nodes from targets that were both sources and targets' % (len(set(sources).intersection(set(targets))))
        targets = [t for t in targets if t not in all_sources]

    print '%d sources and %d targets; %d high penalty sources (empty backward stars)'  % (len(sources),len(targets),len(high_penalty_sources))
    name_dict = {}
    return sources,targets,high_penalty_sources,name_dict


def add_super_nodes(H,sources,targets,high_penalty_sources,dataset):

    '''
    Takes in the sources and targets and creates a supersource and sink
    Returns the supersource and supertarget (sink)
    '''
    super_source = 'SUPERSOURCE'
    H.add_hyperedge(set([super_source]),set(sources.union(high_penalty_sources)),weight=0)
    super_target = 'SUPERTARGET'
    if dataset == 'test' or dataset == 'reactome' or dataset == 'WNT':
        H.add_hyperedge(set(targets),set([super_target]), weight=1)
    else: # dataset == 'ncipid'
        for t in targets:
            H.add_hyperedge(set([t]),set([super_target]), weight=1)

    return super_source,super_target


if __name__ == '__main__':
    main(sys.argv)
