# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 10:38:01 2017

@author: matteo
"""

import sys
from collections import defaultdict

"""
Use case:

grossflow=compute_grossflow(G)

find_persistent_transients MUST be called AFTER find_transients (it needs info on 
    ephemeral transients)
    
    NO - now it computes ephemereal transients by itself (calling eph function)

transients_val_net_grossflow=transients_val/grossflow


In this version, "persistent" split/merge are also looked for; they constitute,
as opposed to "ephemeral" transients, the counterbalancing "interesting" part
to "noise" ("skifo")
"""

def main():
    print "Transients computation library."
    sys.exit(0)

#total nodes size
def compute_grossflow(_G, node_weight_var='size'):

    totSize = sum([ d[node_weight_var] for u,d in _G.nodes(data=True) ])

    return totSize

"""
Ephemereal transients (ETs), <>- and X-kind, are looked for first.
Next, potential persistent transients (splits and merges) are compared to ETs' critical nodes.
If no match occurs, we definitely call them persistent transients (PTs).

Absolute values are computed for both ETs and PTs as sum of the sizes of nodes of interest.
The total sum of the network nodes sizes is used to compute ETs and PTs overall ratios.
"""


def find_pers_ephem_transients(_G, persistence_thr = 0.0, node_weight_var='size', edge_weight_var='weight'):

    m_color_label=set()
    s_color_label=set()
    ms_color_label=set()
    
    print "Transients detection. Persistence threshold =", persistence_thr
    #G=_G.copy()
    G=_G.to_undirected()
    y_range_unsorted=set()

    #find relevant years ("compact" set)
    for u, d in G.nodes(data=True):
        y_range_unsorted.add(d['Xpos'])
        
    y_range=sorted(y_range_unsorted)

    #move forward, look back and forth
    #nodes_to_watch={} useless
    transients={} #both good and bad, nodes: good/bad flag (0 bad, 1 good)
    bad_x_transients_full = {} #central nodes for BAD X transients
    bad_d_transients_full = {} #central node (1) for BAD <> transients
    
    #weird X transients with multiple crossing t2 links, get rid of extra t2 links
    t2_edges_to_go=set()        

    for y in y_range[1:-1]: #we're looking at t +/- 2, censor beginning and ending
        y0_nodes=[u for u, d in G.nodes(data=True) if d['Xpos']==y-1] #previous
        y1_nodes=[u for u, d in G.nodes(data=True) if d['Xpos']==y]   #this year
        y2_nodes=[u for u, d in G.nodes(data=True) if d['Xpos']==y+1]   #next year
        
        successors_y0   = {}
        successors_y1   = {}
        predecessors_y1 = {}
        predecessors_y2 = {}
        #AFAF check unusual (>2 pred/succ) situations, WARNING, there be lions
        for u in y0_nodes:
            _succ_y0  = [n for n in G.neighbors(u) if n in y1_nodes] #of course, only t+1's
            if len(_succ_y0)>=2:
                successors_y0[u]=tuple(sorted(_succ_y0)) #sorted to avoid [x,y] not matching [y,x]

        for u in y1_nodes:
            _succ_y1  = [n for n in G.neighbors(u) if n in y2_nodes] 
            if len(_succ_y1)>=2:
                successors_y1[u]=tuple(sorted(_succ_y1)) #sorted to avoid [x,y] not matching [y,x]

        for u in y1_nodes:
            _pred_y1  = [n for n in G.neighbors(u) if n in y0_nodes] 
            if len(_pred_y1)>=2:
                predecessors_y1[u]=tuple(sorted(_pred_y1))

        for u in y2_nodes:
            _pred_y2  = [n for n in G.neighbors(u) if n in y1_nodes] 
            if len(_pred_y2)>=2:
                predecessors_y2[u]=tuple(sorted(_pred_y2))

        #3 rules for X:
        #i., ii. find couples of y0 AND y2 which share a common y1 node
        #iii. check if there are "t+2" links for ALL pred-succ pairs
        #   more precisely: if we have just 2 pred+2 succ, having BOTH pairs linked by t+2 gves a BAD transient (should not be there)
        #                   if we have MORE than 2 (on either side), we assume it to be bad IIF we have AT LEAST two t+2 links. Deal.

        #i.+ii.: nodes which have both predecessors AND successors: 
        trans_x_centers = [u for u in successors_y1.keys() if u in predecessors_y1.keys()]
        #i.: nodes with multiple successors ONLY
        trans_x_multsucc =  [u for u in successors_y1.keys() if u not in predecessors_y1.keys()]
        #ii.: nodes with multiple predecessors ONLY
        trans_x_multpred = [u for u in predecessors_y1.keys() if u not in successors_y1.keys()]
        
        #first, catch real transients (with both pred and successors)
        #  tell apart + and - transients by checking rule iii.
        trans_x = [ (u,predecessors_y1[u],successors_y1[u]) for u in trans_x_centers ]
        for tr in trans_x:
            focal_node = tr[0]
            _predecessors = tr[1]
            _successors = tr[2]
            tr_w = G.node[focal_node][node_weight_var] #weight (if good)
            transients[focal_node] = tr_w #good            
            print "Making node", focal_node, "good x trans:", transients[focal_node]
            ps_edges = [] #store them here, will need them to re-split central nodes

            for p in _predecessors:     #predecessors (left side)
                for s in _successors: #successors (right side)
                    #print "Checking edge",p,s,"(",tr,")"
                    if G.has_edge(p,s):
                        p_s_edge = G[p][s]
                        #print "... edge",p,s,"(",tr,") is there"
                        
                        p_s_edge_weight = p_s_edge[edge_weight_var]
                        if p_s_edge_weight > persistence_thr:
                            ps_edges.append((p,s,p_s_edge_weight))
                            
            if len(ps_edges) > 1: #see above, 2 or more t+2 links mean bad things
            
                transients[focal_node] *= -1.0 #invert sign
                print "Turning node", focal_node, "bad x trans:", transients[focal_node]

                #also, add to bad X transients list (to later merge)
                bad_x_transients_full[(focal_node,tuple(ps_edges))] = transients[focal_node]

                #WAIT, WAIT: if multiple edges share the same pred/successor, then it's not bad.
                print "Check weirdness in ephemereal lib."
                #idea: when there is a multiple t+2 pre/succ link, keep only the strongest
                e_succ = defaultdict(list)
                e_pred = defaultdict(list)
                for e in ps_edges:
                    u1 = e[0]
                    u2 = e[1]
                    w = e[2]
                    e_succ[u1].append((u2, w))
                    e_pred[u2].append((u1, w))
                    
                for uu,ss in e_succ.items():
                    if len(ss)>1:
                        sorted_succ = sorted(ss, key=lambda x: x[1], reverse = True)
                        #main_e_succ = sorted_succ[0]
                        succ_to_go = sorted_succ[1:]
                        #print "weirdness 1succ_to_go",sorted_succ
                        #print "weirdness 2succ_to_go",succ_to_go
                        for stg,ww in succ_to_go:                    
                            print "popopo", uu,stg
                            t2_edges_to_go.add(tuple(sorted( (uu,stg) )))
                            
                        #print "weirdness 3succ_to_go",t2_edges_to_go
                        
                
                for uu,pp in e_pred.items():
                    if len(pp)>1:
                        sorted_pred = sorted(pp, key=lambda x: x[1], reverse = True)
                        #main_e_pred = sorted_pred[0]
                        pred_to_go = sorted_pred[1:]
                        #print "weirdness 2pred_to_go", pred_to_go
                        for stg,ww in pred_to_go:                            
                            t2_edges_to_go.add(tuple(sorted( (uu,stg) )))
                            
        #print "weirdness 3pred_succ_to_go",t2_edges_to_go

        #3 rules for <>:
        #find candidates (y0 nodes w/multiple successors and y2 w/multiple pred)
        #loop predecessors and successors to find matches (sharing >1 node)
        #i., ii. find couples of y0 AND y2 which share at lest 2 common y1 nodes
        #iii. check if there is a "t+2" link between the pred and the succ.

        #keep track of predecessors, successors implied here, take them out of "leftovers" to use as + transients
        trans_d_implied_nodes_as_splits = []   #merges and splits stubs (originating/ending nodes) must be counted only once, exclude!
        trans_d_implied_nodes_as_merges = []   #merges and splits stubs (originating/ending nodes) must be counted only once, exclude!
        
        for pr, pr_y1s in successors_y0.items():
            for su, su_y1s in predecessors_y2.items():
                intersection = tuple(set(pr_y1s) & set(su_y1s))
                #intersections leading to a one-element set are "Z"'s
                if len(intersection) > 1:
                    #we have a <>
                    trans_d_implied_nodes_as_splits.append(pr)
                    trans_d_implied_nodes_as_merges.append(su)
                    
                    tr_w = sum( [G.node[u][node_weight_var] for u in intersection] )
                    transients[intersection]=tr_w
                    
                    #check t+2 link
                    if G.has_edge(pr,su):
                        p_s_edge = G[pr][su]
                        if p_s_edge[edge_weight_var] > persistence_thr:
                            print "Making nodes", intersection, "bad <> trans because a link between",pr,"and",su,"exists."
                            transients[intersection] *= -1.0 #invert sign
                            #also, add to bad X transients list (to later merge)
                            bad_d_transients_full[(intersection)] = transients[intersection] #already negative
                            
                            #if bad, starting and ending nodes are not good for previous merges or successive splits, either
                            trans_d_implied_nodes_as_merges.append(pr)
                            trans_d_implied_nodes_as_splits.append(su)
                            
                if len(intersection) > 2:
                    print "Warning: <> transient w/ more than 2 shared nodes:", intersection

        
        #nodes with multiple successors only (<) not excluded
        nodes_good_stubs1 = [u for u in (trans_x_multsucc) if u not in trans_d_implied_nodes_as_splits]
        nodes_good_stubs2 = [u for u in (trans_x_multpred) if u not in trans_d_implied_nodes_as_merges]

        nodes_good_stubs1_w_successors = {}        
        for k,v in successors_y1.items():
            if k in nodes_good_stubs1:
                nodes_good_stubs1_w_successors[k]=v

        nodes_good_stubs2_w_predecessors = {}        
        for k,v in predecessors_y1.items():
            if k in nodes_good_stubs2:
                nodes_good_stubs2_w_predecessors[k]=v

        print "NODES GOOD STUBS w/succ", nodes_good_stubs1_w_successors
        
        #nodes_good_stubs = nodes_good_stubs1+nodes_good_stubs2
        #print "NODES GOOD STUBS", nodes_good_stubs
        #nodes_good_stubs2 = [u for u in (trans_x_multpred + trans_x_multsucc) if u not in nodes_in_bad_transients]
        
        #one loop for splits, one for merges (important nodes differ)
        for u, s in nodes_good_stubs1_w_successors.items() : #good splits, highlight successors nodes (will look like one transient each)
            #tr_w_sum = 0.0
            print "COLOR - good split;", u ,s
            for v in s:
                tr_w = G.node[v][node_weight_var] #weight (always good)
                transients[v] = tr_w
                print "COLOR2 - good split;", v
                s_color_label.add(v)

            #tr_w = tr_w_sum

        for u in nodes_good_stubs2 : #good merges, highlight the one successor node
            tr_w = G.node[u][node_weight_var] #weight (always good)
            transients[u] = tr_w
            
        #prov for coloring
        for u, p in nodes_good_stubs2_w_predecessors.items() : #good merges, highlight the one successor node
            print "COLOR - good merge;", u, p
            print "COLOR2 - good merge;", u
            m_color_label.add(u)


    ms_color_label = m_color_label.intersection(s_color_label)
    mo_color_label = m_color_label.difference(ms_color_label)
    so_color_label = s_color_label.difference(ms_color_label)

    """
    for _n in mo_color_label:
        print "COLOR3 MERGE;", idfromlabel(G, _n), "; merge"
    for _n in so_color_label:
        print "COLOR3 SPLIT;", _n, "; split"
    for _n in ms_color_label:
        print "COLOR3 BOTH ;", _n, "; both"
    """
    
    tr_total_weight = sum(transients.values())
    #split negative and positive transients
    pos_trans = [_w for _w in transients.values() if _w>0]
    neg_trans = [_w for _w in transients.values() if _w<0]

    print "TRANSIENTS from ephem lib:", y, sum(transients.values()),"(",sum(pos_trans),"/",sum(neg_trans),")"

       
    return transients, tr_total_weight,\
        bad_x_transients_full, bad_d_transients_full,\
        t2_edges_to_go 


def idfromlabel(_G, lab):
    #assume label uniqueness
    _id = [k['id'] for u, k in _G.nodes(data=True) if k['label'] == lab]

    return _id
    
if __name__ == "__main__":
    main()
