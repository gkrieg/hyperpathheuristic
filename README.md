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

Currently under construction...


