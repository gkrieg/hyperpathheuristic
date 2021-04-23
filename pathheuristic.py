import itertools
import heapq

from halp.directed_hypergraph import DirectedHypergraph
from halp.algorithms.directed_paths import b_visit

def tail_path_heuristic(H,source,target,node_dict={}):
    '''
    Finds the distance to the tail for each edge in the hypergraph H that is reachable from the source vertex and backwards-recoverable from the target
    Returns a list of the tail distance for each edge
    ''' 
    reachableedges,edgedict,taildistancelist,heap,reachableedgecounter,reachedtable,entry_finder,counter,H = initialize(H,source,target,node_dict=node_dict)

    returnpath = []
    while reachableedgecounter > 0:
        e,epriority = pop_node(heap,entry_finder)
        
        if e in reachableedges:
            reachableedgecounter -= 1
            reachableedges.remove(e) #possibly don't need this
            edgedict[e]['isremoved'] = True
            a,b,c = recover(H,e,'full',edgedict)
            P,inedges = trim(a,b,c,source,taildistancelist)
            if target in H.get_hyperedge_head(e):
                returnpath = P
                for p in P:
                    tail = []
                    unreadabletail = []
                    for v in H.get_hyperedge_tail(p):
                        if len(node_dict) > 0:
                            tail.append(node_dict[v])
                            unreadabletail.append(v)
                        else:
                            tail.append(H.get_node_attribute(v,'label'))
                            unreadabletail.append(v)
                    head = []
                    unreadablehead = []
                    for v in H.get_hyperedge_head(p):
                        if len(node_dict) > 0:
                            head.append(node_dict[v])
                            unreadablehead.append(v)
                        else:
                            head.append(H.get_node_attribute(v,'label'))
                            unreadablehead.append(v)

            taildistancelist[e] = weight(H,P) - weight(H,[e])
            edgedict[e]['bestinedges'] = inedges
            edgestoprocess = findedgestoprocess(e,H,reachedtable,edgedict)
            for f in edgestoprocess:
                edgedict[f]['candidateinedges'].append(e)
                if f in entry_finder and edgedict[f]['isremoved'] == False:
                    a,b,c = recover(H,f,'short',edgedict)
                    fpath,_ = trim(a,b,c,source,taildistancelist)
                    edgepriority = weight(H,fpath)
                    add_node(heap,f,edgepriority,counter,entry_finder)
                elif edgedict[f]['tailcount'] == 0 and f not in entry_finder:
                    a,b,c = recover(H,f,'short',edgedict)
                    fpath,_ = trim(a,b,c,source,taildistancelist)
                    edgepriority = weight(H,fpath)
                    add_node(heap,f,edgepriority,counter,entry_finder)
        else:
            print('edge {} was not in the reachable set but still in the hypergraph'.format(e))

    return taildistancelist,H,returnpath


def findedgestoprocess(e,H,reachedtable,edgedict):
    '''
    Finds the edges whose tails intersect with e's head. Also sets e to their inedge lists and decrements their reached counter in the reachedtable
    returns the list of edges to process
    '''
    F = set()
    for v in H.get_hyperedge_head(e):
        for f in H.get_forward_star(v):
            if v not in reachedtable[f]:
                reachedtable[f].append(v)
                edgedict[f]['tailcount'] -= 1
                F.add(f)
    return list(F)


def weight(H,F):
    '''
    Takes as input a set of edges and returns the sum of the weight of the edges
    returns the weight of the edge set
    '''
    w = 0
    for f in F:
        w += H.get_hyperedge_weight(f)
    return w

def recover(H,e,flag,edgedict):
    '''
    finds the set of edges that are recovered recursively for each edge starting with the inedge list of e. Returns this set of edges
    The flag sets whether to use bestinedges or candidate inedges for every edge except e, or to use all inedges for all edges
    Returns the set of edges backward recoverable from the edge
    '''
    Q = []
    Qindex = 0
    F = set()
    if flag == 'all':
        for v in H.get_hyperedge_tail(e):
            for f in H.get_backward_star(v):
                if f not in F:
                    Q.append(f)
                    F.add(f)
    if flag != 'all':
        for f in edgedict[e]['candidateinedges']:
            Q.append(f)
            F.add(f)
    while Qindex < len(Q):
        #process the next edge in the Q
        f = Q[Qindex]
        Qindex += 1
        if flag == 'all':
            for v in H.get_hyperedge_tail(f):
                for g in H.get_backward_star(v):
                    if g not in F:
                        F.add(g)
                        Q.append(g)
        elif flag == 'full':
            for g in edgedict[f]['candidateinedges']:
                if g not in F:
                    F.add(g)
                    Q.append(g)
        elif flag == 'short':
            for g in edgedict[f]['bestinedges']:
                if g not in F:
                    F.add(g)
                    Q.append(g)
        else:
            print('something wrong in recover with the flag')

    return H,F,e

def trim(H,F,e,source,taildistancelist):
    '''
    Takes the edge list F, which is a super path from the source vertex to e and repeatedly removes an edge and checks reachability for each edge in F
    Returns the trimmed edge list and prints that the sink is not reachable from the source if the original edge set is not a superpath
    '''
    sortededges = []
    if e[0] != 'e':
        reachingset = [e]
    else:
        reachingset = H.get_hyperedge_tail(e)
        
    for f in F:
        sortededges.append((taildistancelist[f],f))
    sortededges.sort()
    justedges = [s[1] for s in sortededges]
    H2 = DirectedHypergraph()
    H2.add_node(source)
    H2edgeids = {}
    for f in justedges:
        newf = H2.add_hyperedge(H.get_hyperedge_tail(f),H.get_hyperedge_head(f))
        H2edgeids[f] = newf
    nodeset,_,__,___ = b_visit(H2,source)
    isereached = True
    for v in reachingset:
        if v not in nodeset:
            isereached = False
    if isereached == False:
        print('invalid edge set given to trim. Sink not reachable from source')
        print(nodeset)
    F2 = []
    tailedges = []
    for f in justedges:
        H2.remove_hyperedge(H2edgeids[f])
        nodeset,_,__,___ = b_visit(H2,source)
        isereached = True
        for v in reachingset:
            if v not in nodeset:
                isereached = False
        if isereached == False:
            #This means we cannot remove f
            newf = H2.add_hyperedge(H.get_hyperedge_tail(f),H.get_hyperedge_head(f))
            H2edgeids[f] = newf
            F2.append(f)
            if e[0] == 'e':
                for v in H.get_hyperedge_tail(e):
                    if v in H.get_hyperedge_head(f):
                        tailedges.append(f)
    if e[0] == 'e':
        F2.append(e)
    return F2,tailedges

def initialize(H,source,target,node_dict={}):
    '''
    Finds the set of reachable and backwards-recoverable edges, sets up the tail reached counters, the inedge lists,  and the heap pointer for each edge.
    Initializes the taillength for each edge, and the heap, and the reachable edge counter
    '''
    edgedict = {}

    #find reachable and backwards-recoverable edge list
    reachableedges = findreachableandbackrecoverable(H,source,target)

    #trim H
    H2 = DirectedHypergraph()
    if len(node_dict) == 0:
        for v in H.node_iterator():
            H2.add_node(v,H.get_node_attributes(v))
    for edge in reachableedges:
        H2.add_hyperedge(H.get_hyperedge_tail(edge),H.get_hyperedge_head(edge),weight=H.get_hyperedge_weight(edge))
    H2.add_node('SUPERSOURCE',{'label': 'SUPERSOURCE'}) 
    for edge in H2.get_forward_star('SUPERSOURCE'):
        forwardstarlist = []
        for v in H2.get_hyperedge_head(edge):
            if len(H2.get_forward_star(v)) > 0:
                forwardstarlist.append(v)
        H2.remove_hyperedge(edge)
        if len(forwardstarlist) > 0:
            H2.add_hyperedge(['SUPERSOURCE'],forwardstarlist,weight=0)
        else:
            H2.add_hyperedge(['SUPERSOURCE'],[],weight=0)

    H = H2
    for v in H.get_node_set():
        #need to remap the edges because just calling remove_node(v) removes also all the hyperedges v is in
        if len(H.get_forward_star(v)) == 0 and v != 'SUPERTARGET':
            #find edges that need to be replaced
            backedges = H.get_backward_star(v)
            for e in backedges:
                tail = H.get_hyperedge_tail(e)
                head = H.get_hyperedge_head(e)
                head.remove(v)
                w = H.get_hyperedge_weight(e)
                H.remove_hyperedge(e)
                H.add_hyperedge(tail,head,weight=w)
            H.remove_node(v)
    reachableedges = []
    H2.add_node('SUPERSOURCE',{'label': 'SUPERSOURCE'}) 
    for edge in H.hyperedge_id_iterator():
        reachableedges.append(edge)

    #initialize edgedict
    for e in H.hyperedge_id_iterator():
        edgedict[e] = {'isremoved': False, 'bestinedges': [], 'candidateinedges': [], 'tailcount': len(H.get_hyperedge_tail(e))}

    #initialize taildistancelist
    taildistancelist = {}
    for e in reachableedges:
        taildistancelist[e] = 'inf'

    #initialize reachableedgecounter
    reachableedgecounter = len(reachableedges)

    #initialize heap
    heap = []
    entry_finder = {}               # mapping of tasks to entries, this and the line below are strictly for the heap
    counter = itertools.count()     # unique sequence count
    for e in H.get_forward_star(source):
        if e in reachableedges:
            add_node(heap,e,weight(H,[e]),counter,entry_finder)

    #initialize reachedtable
    reachedtable = {}
    for e in H.hyperedge_id_iterator():
        reachedtable[e] = []

    return reachableedges,edgedict,taildistancelist,heap,reachableedgecounter,reachedtable, entry_finder, counter, H

def findreachableandbackrecoverable(H,source,target):
    '''
    Finds and returns the set of reachable and backwards-recoverable edges
    '''
    reachable = findreachable(H,source)

    recoverable = findrecoverable(H,target)
    
    reachrecover = list(set(reachable) & set(recoverable))
    return reachrecover

def findreachable(H,source):
    '''
    Finds the set of reachable edges from the source
    '''
    _,__,edges,___ = b_visit(H,source)
    reachable = []
    for e in edges:
        if edges[e] != None:
            reachable.append(e)
    
    return reachable

def findrecoverable(H,target):
    '''
    Finds the set of backwards-recoverable edges from the sink
    '''
    sinkedges = []
    for e in H.get_backward_star(target):
        sinkedges.append(e)

    recoverable = set(sinkedges)
    for e in sinkedges:
        _,R,__ = recover(H,e,'all',{})
        for f in R:
            recoverable.add(f)
    return recoverable
        

###################################################################################
                                #heapq helper functions
###################################################################################

def add_node(pq,node, priority, counter, entry_finder):
    'Add a new task or update the priority of an existing task'
    if node in entry_finder:
        remove_node(pq,node,entry_finder)
    count = next(counter)
    entry = [priority, count, node]
    entry_finder[node] = entry
    heapq.heappush(pq, entry)

def remove_node(pq,node,entry_finder):
    'Mark an existing task as REMOVED.  Raise KeyError if not found.'
    REMOVED = '<removed-node>'      # placeholder for a removed node
    entry = entry_finder.pop(node)
    entry[-1] = REMOVED

def pop_node(pq,entry_finder):
    'Remove and return the lowest priority task. Raise KeyError if empty.'
    REMOVED = '<removed-node>'      # placeholder for a removed node
    while pq:
        priority, count, node = heapq.heappop(pq)
        if node != REMOVED:
            del entry_finder[node]
            return node,priority
    raise KeyError('pop from an empty priority queue')

def emptypq(pq,entry_finder):
    REMOVED = '<removed-node>'      # placeholder for a removed node
    while pq:
        priority, count, node = heapq.heappop(pq)
        if node != REMOVED:
            heapq.heappush(pq,[priority, count, node])
            return False
    return True
