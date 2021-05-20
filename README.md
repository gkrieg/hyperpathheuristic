# hyperpaths
## Prerequisites

The code for constructing directed hypergraph objects (required for the hyperpath heuristic, cut finding procedure, and cutting planes algorithm) is available at https://github.com/annaritz/pathway-connectivity.
Once you have converted the OWL file from Reactome or Pathway Commons, you should have a -hypernodes.txt file and a -hyperedges.txt file, which will be read in by run.py to construct the directed hyperpath object.

To run the hyperpath heuristic and hyperpath generation algorithm, only python is required, as well as halp, which can be installed via pip.

For the cutting-planes algorithm, the LP solver CPLEX is required, which can be installed through https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html

## Installation

To install, you must create a folder in the directory called parsed/ and then set ROOTDIR in run.py to be the current directory. It is in this parsed/ directory where your hypernodes.txt and hyperedges.txt files should go.

## Algorithms

The hyperpath heuristic can be run using the following command:

```
 python run.py --heuristic --source <source> --target <target> --name <prefix of hypernodes.txt and hyperedges.txt file>
```

The hyperpath enumeration algorithm can be run with the following command

```
 python run.py --enumeration --source <source> --target <target> --name <prefix of hypernodes.txt and hyperedges.txt file> 
```

## Datasets

Parsed versions of the hypergraphs built off of NCI-PID and Reactome datasets are available in parsedhypergraphs. These hypergraphs can be searched from run.py. We built hypergraphs using three subsets of NCI-PID pathways, and call these the Small (WNT5A), Medium (WNT), and Large(allpid) datasets. The Small dataset is a small Wnt signaling pathway consisting of the union of two pathways: “degradation of β-catenin” and “canonical Wnt signaling”. The Medium dataset is a larger Wnt signaling pathway including four additional pathways: “noncanonical Wnt signaling”, “Wnt signaling network”, “regulation of nuclear β-catenin”, and “presenilin action in Notch and Wnt signaling”, which correspond to non-canonical branches of Wnt signaling. The Large dataset contains all NCI-PID pathways. Similarly, the Reactome dataset is the union of all Reactome pathways. The NCI-PID and Reactome pathways were downloaded in the BioPAX format from Pathway Commons and processed using a parser referenced above.

## Results

The results directory contains two files which summarize the instances from the datasets where our heuristic returns a cyclic hyperpath. Each file contains instances in the form of:

<pre>
line with target id
line with cyclic vertex id
description of path:
    one line for each edge, giving the screen name from Pathway Commons of each vertex in the tail then head of the hyperege
</pre>

Currently under construction...


