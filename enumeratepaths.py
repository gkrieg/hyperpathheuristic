import pickle as pkl
from halp.directed_hypergraph import DirectedHypergraph
from pathheuristic import weight,tail_path_heuristic

from collections import deque
def enumerateallpaths(H,P,source,target,node_dict=None):
    Q = deque()
    #represent problems by the edges they exclude from the hypergraph H
    inlist = []
    for e in P:
        if 'SUPERSOURCE' not in H.get_hyperedge_tail(e) and 'SUPERTARGET' not in H.get_hyperedge_head(e):
            Q.append((set(inlist),set([e])))
            inlist.append(e)
            print(e)

    #now that the queue is seeded, repeatedly pop, solve, append with new edges
    #If no path exists, this means we do not append new problems
    #We will also keep it so that you will only ever add edges to the set with a higher edge label, never lower
    bestweight = weight(H,P)
    bestP = P
    paths = set([tuple(sorted(P))])

    while len(Q) > 0:
        inlist,removededges = Q.popleft()
        H2 = H.copy()
        for e in removededges:
            H2.remove_hyperedge(e)

        _,H3,Ptemp = tail_path_heuristic(H2,source,target,node_dict=node_dict)

        if Ptemp != []:
            w = weight(H3,Ptemp)
            if w < bestweight:
                bestweight = w
                bestP = Ptemp
            Hpath = []
            for e in Ptemp:
                # for each edgein order, you put one in the out list you inherited from this problem and inherit the in list from the parent. The first subproblem removes e1. The second puts e1 in the in list and e2 in the outlist, etc.
                #If the edge is in your in list it cannot be added to your outlist
                if 'SUPERSOURCE' in H3.get_hyperedge_tail(e) or 'SUPERTARGET' in H3.get_hyperedge_head(e):
                    continue
                newe = H.get_hyperedge_id(H3.get_hyperedge_tail(e),H3.get_hyperedge_head(e))
                Hpath.append(newe)
                if newe not in inlist:
                    newremovededges = removededges.union(set([newe]))
                    print('new removed edges',newremovededges)
                    print(newe)
                    Q.append((inlist,newremovededges))
                    inlist.add(newe)
            paths.add(tuple(sorted(Hpath)))
    return paths,bestP,bestweight


def num(e):
    return int(e[1:])





