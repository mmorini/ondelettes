import networkx as nx
from operator import itemgetter
import sys
import os
import scipy.stats as ss
from random import random #shuffle, random
from copy import deepcopy #for nodes attributes, otherwise... bang

from networkx.readwrite import json_graph
import json

import community
import ephemereal_3rules
import histUtils as hu
import IDCards_maker as idcm

inpath='/home/matteo/WORK/PhD/projets/ONDELETTES/FINAL/4/scripts'
outpath='/home/matteo/WORK/PhD/projets/ONDELETTES/FINAL/4/scripts'

filename = 'g_big_w.gexf'


firstNodes = 10   #filt 1: abs number per year
percNodes = .1    #filt 2: % of nodes per year
minSize = 10      #filt 3: min node size
minCommPerc = .01 #filt 4: min community percentage on yearly total

FILTER = 4 ###################################

def dump_net_json(G, myfile):
    _data = json_graph.node_link_data(G)
    #json_data = json.dumps(_data)
    json_filename = myfile+'.json'
    with open(json_filename, 'w') as outfile:
        json.dump(_data, outfile)
    return

def split_node_name(u):
    #print "U",u
    splt = u.split('_')
    year = int(splt[1]) #ugly hack, but who cares
    comm = int(splt[0]) 
    
    return year, str(comm)

def change_comm_name_merge(cn):
    if len(cn)<6:
        ncn = cn+'00000'
    else:
        ncn = cn
    return ncn

def change_comm_name_split(cn, n):
    if len(cn)<6:
        if cn=='0': #fix for BC, cannot deal with zero-leading communities
            cn='999'
        ncn = cn+'0000'+str(n)
    else:
        ncn = cn
    return ncn

def part2nodes(part):
    u_dict = {}
    for p, us in part.items():
        for u in us:
            u_dict[u]=p
            
    return u_dict

"""
def avgMaxJac(_G):
    #find first/last year
    y1 = min([k['year'] for u,k in _G.nodes(data=True)])
    y2 = max([k['year'] for u,k in _G.nodes(data=True)])

    print y1,y2

    unique_edges=set()
    unique_nodes=set()
    avgmaxjac_tot = 0.0
    
    for y in range(y1+1,y2): #no last year, ok
        y0_nodes = [u for u,k in _G.nodes(data=True) if k['year']==y-1]
        y1_nodes = [u for u,k in _G.nodes(data=True) if k['year']==y]
        y2_nodes = [u for u,k in _G.nodes(data=True) if k['year']==y+1]

        avgmaxjac_num = 0.0
        avgmaxjac_den = sum( [ k['size'] for u,k in _G.nodes(data=True) if u in y1_nodes])    
        
        #for every node keep only the strongest link
        for u in y1_nodes:
            _nei_edges_0 = [(u, v) for v in _G.neighbors(u) if v in (y0_nodes)]
            _nei_edges_2 = [(u, v) for v in _G.neighbors(u) if v in (y2_nodes)]

            _nei_edges_wdata_0 = [ (u, v, _G.edge[u][v]['weight']) for u,v in _nei_edges_0 ]
            _nei_edges_wdata_2 = [ (u, v, _G.edge[u][v]['weight']) for u,v in _nei_edges_2 ]
            
            if len(_nei_edges_wdata_0)>0:
                max_e_0 = max( _nei_edges_wdata_0, key = itemgetter(2))
                _e = (max_e_0[0], max_e_0[1])
                if _e not in unique_edges: #count every edge only once
                    unique_edges.add(_e)
                    avg_comm_size = (_G.node[max_e_0[0]]['size']*1.0) + (_G.node[max_e_0[1]]['size']*1.0)/2.0
                    max_jac_0 = avg_comm_size * max_e_0[2] #edge weight == jaccard
                    avgmaxjac_num += max_jac_0

            if len(_nei_edges_wdata_2)>0:
                max_e_2 = max( _nei_edges_wdata_2, key = itemgetter(2))
                _e = (max_e_2[0], max_e_2[1])
                if _e not in unique_edges: #count every edge only once
                    unique_edges.add(_e)
                    avg_comm_size = (_G.node[max_e_2[0]]['size']*1.0) + (_G.node[max_e_2[1]]['size']*1.0)/2.0
                    max_jac_2 = avg_comm_size * max_e_2[2]           
                    avgmaxjac_num += max_jac_2
                    
        avgmaxjac_year = avgmaxjac_num*1.0 / avgmaxjac_den*1.0
        print avgmaxjac_year
        avgmaxjac_tot += avgmaxjac_year
    

    return avgmaxjac_tot
"""

SCHIFO_BY_ALPHA={}

myfile=(os.path.join(inpath,filename))

g = nx.read_gexf(myfile)

#set t2 links aside, bring them back after the process
t2_links = [(u,v,d) for u,v,d in g.edges(data=True) if d['scope']=='t2']

print "T2 links:", t2_links

for edge in g.edges():
    if len(g.get_edge_data(edge[0], edge[1])) == 1: #no weight??
        g.get_edge_data(edge[0], edge[1])['weight'] = 1

#years = range(min( [g.node[node]['year'] for node in g] ), max( [g.node[node]['year'] for node in g] ))# + 1) #list of years
years = range(min( [g.node[node]['year'] for node in g] ), max( [g.node[node]['year'] for node in g] ) + 1) #list of years
print "YEARS:",years

### From right to left
### let's check mutual strongest edge
### Weights are NOT normalized
### Define all the trends
### keep track of all edges between different trends
### rearange from the lowest to the highest edges


C = {}
C_static = {}
Y = {}
for year in years:
    C[year] = [node for node in g if g.node[node]['year'] == year]
    C_static[year] = [node for node in g if g.node[node]['year'] == year]
    Y[year] = []   
#dictionary C contains for every year list of communities that belong to that year
#community with y-pos assigned will be removed from the list

#dictionary Y is now empty but will contain taken y-positions for every year

brojac = 0 # brojac counts the order of the trend
trends = {} # will contain all the trends

def move_down_from_trend(com):
    for node in g:
        if g.node[node]['y-pos'] < g.node[com]['y-pos']:
            g.node[node]['y-pos'] = g.node[node]['y-pos'] - 1
sklj = 0
def reposition_trends(trend_big, trend_small):
    com = trend_big[0]
    move_down_from_trend(com)
    for node in trend_small:
        g.node[node]['y-pos'] = g.node[com]['y-pos'] - 1
        g.node[node]['came_from'] = com
        g.node[node]['order_of_position'] = sklj

#Y[max(years)].append(len(C[max(years)]) * 10) #first y position, for the last year's community

while True:
    ### starting a trend
    trend = [] #this list will contain a current trend, just temporaly
    brojac += 1
    current_year = max([year for year in years if C[year] != []]) #start from the last year that contains unassigned communities
    tmp = max([g.node[node]['size'] for node in C[current_year]])
    C1 = [node for node in C[current_year] if g.node[node]['size'] == tmp]
    C1 = C1[0] #not assigned community with maximum size for the last year
    ### with this C1 the trend starts  
    
    g.node[C1]['y-pos'] = brojac #this is just temporal position, every trend has their own y position
    C[current_year].remove(C1) #remove community with assigned y-position from the list of the communities for the current_year
    trend.append(C1)
    
    #while current_year > min(years) :
    while ( current_year > min(years) and current_year < max(years)):
        
        ### here we make the trend that begins with C1
        
        current_year -= 1 #now look in the previous year
        intsect_left = list(set(g.neighbors(C1)) & set(C[current_year])) #list of communities from current year that have edges with C1 from next year
        if intsect_left == []:
            break
        tmp = max([g.get_edge_data(C1, node)['weight'] for node in intsect_left])
        C2 = [node for node in intsect_left if g.get_edge_data(C1, node)['weight'] == tmp]
        C2 = C2[0] #community with maximum edge weight with C1
        
        #print "---", C2
        intsect_right = list(set(g.neighbors(C2)) & set(C_static[current_year + 1])) #chosen community C2, check right neighbors weights
        tmp2 = max([g.get_edge_data(C2, node)['weight'] for node in intsect_right]) #maximum weight for C2
        
        if tmp2 > tmp: #if C1 C2 is not mutually highest weight
            break
            
        g.node[C2]['y-pos'] = g.node[C1]['y-pos'] #assign to this community the same y-pos as C1
        trend.append(C2)
        C[current_year].remove(C2) #remove assigned community
        #Y[current_year].append(g.node[C2]['y-pos']) #add assigned y-position to the list
        C1 = C2

    trends[brojac] = trend    
    if sum(len(C[year]) for year in years) == 0:
        break
    
### at this point trends are made and kept in dict trends
# nx.write_gexf(g, graph_path + "g_w_position_test0.gexf")

### now take the information about heighest weights between trends

trends_edges = []

#print "TRENDS", trends.items()

for i in trends.keys(): #go through all the trends
    for node in trends[i]:
        # checking the connections to other communities but not to the one in this trend
        for node2 in list(set(g.neighbors(node)) - set(trends[i])):
            #print "++++++++++++"
            #print "1", g.node[node]
            #print "2", g.node[node2]
            if g.node[node]['y-pos'] < g.node[node2]['y-pos']:
                trends_edges.append([g.node[node]['y-pos'], g.node[node2]['y-pos'], g.get_edge_data(node, node2)['weight']])
        
### now arange trends, using weights between them, starting with lowest weights towards highest
trends_edges.sort(key = itemgetter(2))

for triple in trends_edges: #trends[triple[0]], trends[triple[1]] with weight triple[2]
    if len(trends[triple[0]]) > len(trends[triple[1]]):
        sklj += 1
        reposition_trends(trends[triple[0]], trends[triple[1]])
    else:
        sklj += 1
        reposition_trends(trends[triple[1]], trends[triple[0]])  

pozicije = []
for node in g:
    pozicije.append(g.node[node]['y-pos'])
pozicije = list(set(pozicije))
pozicije.sort()

for i in range(len(pozicije)):
    for node in [node for node in g if g.node[node]['y-pos'] == pozicije[i]]:
        g.node[node]['y-position'] = i
for node in g:
    g.node[node]['y-pos'] = g.node[node]['y-position']
tmp = min([g.node[node]['y-pos'] for node in g])
for node in g:
    g.node[node]['y-pos'] = g.node[node]['y-pos'] - tmp
    del g.node[node]['y-position']


#compute and sort streams by lenghts (provisionally, then jaccard cross-matching)
trend_weights={}
for t, nlist in trends.items():
    trend_weights[t]=sum([d['size'] for u, d in g.nodes(data=True) if u in trends[t]])
sorted_trends_w_weights = sorted(trend_weights.items(), key=itemgetter(1), reverse=True) #by value, aka size
sorted_trends = [t[0] for t in sorted_trends_w_weights]

#reposition nodes according to trends
reverse_trend_dict={}
for t, nlist in trends.items():
    for n in nlist:
        reverse_trend_dict[n]=t
    
for node in g:
    #position at the n-th ranked trend by weight
    g.node[node]['y-pos'] = sorted_trends.index(reverse_trend_dict[node])

#keep max NY nodes per year    , prune edges which are not predec/success
g_filt=nx.Graph()


#filter communities (nodes) according to chosen rule
nodesToKeep=[]
for year in years:
    #list of nodes and sizes for given year
    nodes_year = {node: g.node[node]['size'] for node in g if g.node[node]['year'] == year}
    num_nodes_year = len(nodes_year)
    totNodeSize_year = sum(nodes_year.values())
    #print year, totNodeSize_year

    if FILTER == 1:        
        first_n_nodes_w_size = sorted(nodes_year.items(), key=itemgetter(1), reverse=True)[0:firstNodes] #abs by value, aka size
        first_n_nodes = [n[0] for n in first_n_nodes_w_size]
    elif FILTER == 2:
        num_nodes_to_keep = int(num_nodes_year*percNodes)
        first_n_nodes_w_size = sorted(nodes_year.items(), key=itemgetter(1), reverse=True)[0:num_nodes_to_keep] #% by val, aka size
        first_n_nodes = [n[0] for n in first_n_nodes_w_size]
    elif FILTER == 3:
        first_n_nodes = [u for u, d in nodes_year.items() if d >= minSize ]
    elif FILTER == 4:
        minCommSize = totNodeSize_year*minCommPerc
        first_n_nodes = [u for u, d in nodes_year.items() if d >=  minCommSize]
    else:
        print "INVALID FILTERING RULE."
        sys.exit()

    print "Year %s: %s communities left (out of %s) after filtering, %f %%." %\
        ( year, len(first_n_nodes), len(nodes_year.keys()), (len(first_n_nodes)*100.0/len(nodes_year.keys())) )
    
    nodesToKeep.extend(first_n_nodes)
    
#reposition vertically new subgraph
g_filt = nx.Graph(nx.subgraph(g, nodesToKeep))
nodes = g_filt.nodes()
nodes_data=g_filt.nodes(data=True)

ypositions=[d['y-pos'] for u,d in nodes_data]

ypositions_ranks=ss.rankdata(ypositions, method='dense')

ypos_rank_list=list(ypositions_ranks)

node_yposfix={}
for u in nodes:
    node_yposfix[u] = int(ypos_rank_list.pop(0))

for u,d in g_filt.nodes(data=True):
    g_filt.node[u]['y-pos-fix']=node_yposfix[u]
for t,nlist in trends.items():
    for n in nlist:
        if n in nodes:
            g_filt.node[n]['color']=t


nx.write_gexf(g_filt, myfile+'_layedout.gexf')

    
#add viz position, color, size specifications

g = nx.read_gexf(myfile+'_layedout.gexf')

#bring back t2 edges

t1_nodes=g.nodes()
t2_links_onlyexisting = [(u,v,d) for (u,v,d) in t2_links if (u in t1_nodes and v in t1_nodes)]
g.add_edges_from(t2_links_onlyexisting)

nx.write_gexf(g, myfile+'_viz.gexf')

#also 
dump_net_json(g, myfile)

transients, tot_trans_weight, _, bad_d_transients_full, _ = ephemereal_3rules.find_pers_ephem_transients(g, persistence_thr=0.0)


tot_graph_weight = ephemereal_3rules.compute_grossflow(g)
graph_score = (tot_trans_weight*1.0) / (tot_graph_weight*1.0)
print "ORIGINAL Graph %s score: %f" % (filename, graph_score)

##########################################################  
##########################################################  
### CLEANING UP TRANSIENTS ###############################
##########################################################  
##########################################################  

#make a dict structure, easy to keep updated, with years/communities/nodes
history_master={}

for y in years:
    [commsize_y0, arts_y0, G_y0] = hu.read_gexf_year(y)
    history_master[y] = [commsize_y0, arts_y0, G_y0]
    
####compute scores (skifo, modularity) for PRE GRAPHS
for y in years:
    [_, _part, _G] = history_master[y] 
    _part_as_nodes_comms = part2nodes(_part)
    _mod = community.modularity(_part_as_nodes_comms, _G)
    print "PRE-SCORES", y, _mod

    
    
#print history_master
#AFAF substitute successive calls to hu with history master dictionary

Gcopy=g.copy()
#SORT by descending weight? <<< AFAFAF
#iteratively, fuse all pairs of bad <> transient nodes

#sort by descending magnitude
sorted_bad_d_transients_full = list( sorted( bad_d_transients_full.items(), key=itemgetter(1) ) )

#sort by ascending magnitude
#sorted_bad_d_transients_full = list( sorted( bad_d_transients_full.items(), key=itemgetter(1), reverse = True ) )

#shuffle
#shuffle(sorted_bad_d_transients_full )

#check untreated transients
#untreated_transients=set()
x_weird_transients=set() #not bad transients ("weird"), exclude them


#get rid of all bad <> transients

it=0
prev_trans=[]
for bt, w in sorted_bad_d_transients_full:

    #if len(bt)==2:
    print "Merging %i nodes" % (len(bt))
    merging_node = bt[0]
    
    mn_year, mn_comm = split_node_name(merging_node)  
    mn_size = Gcopy.node[merging_node]['size']

    aggr_size = mn_size
    for on in bt[1:]:
        v = on
        _, v_comm = split_node_name(v)
        
        u_d = '-'.join((Gcopy.node[merging_node]['mainauth'], Gcopy.node[merging_node]['keyword']))
        v_d = '-'.join((Gcopy.node[v]['mainauth'], Gcopy.node[v]['keyword']))

        #merge happens here. Beware of size
        v_size = Gcopy.node[v]['size']
        aggr_size += v_size #accumulate 
        
        Gcopy = nx.contracted_nodes(Gcopy, merging_node, v) #cotracted_node is "u", v disappears
        Gcopy.node[merging_node]['size'] = aggr_size 
        
        #update history master 
        #[0] - {comm: size}
        #[1] - {comm: arts list}
        #[2] - {G}
        #-> nodes merged
        
        #[1] fix resulting merged community
        new_comm_label = change_comm_name_merge(mn_comm)

        #problem: community name has changed with trailing zeros

        if mn_comm not in history_master[mn_year][1]:
            nodes_from_mn = history_master[mn_year][1][change_comm_name_merge(mn_comm)]
        else:
            nodes_from_mn = history_master[mn_year][1][mn_comm]
            del history_master[mn_year][1][mn_comm]
        
        nodes_from_v = history_master[mn_year][1][v_comm]
        nodes_all = nodes_from_mn + nodes_from_v
                
        history_master[mn_year][1][new_comm_label] = nodes_all #add new node
        #[0]        
        history_master[mn_year][0][new_comm_label] = len(nodes_all)

        del history_master[mn_year][1][v_comm]

        #[2]
        #history_master[mn_year][2]=Gcopy
        #print "DDD", mn_year, "DDD", Gcopy.nodes()

        transients, tot_trans_weight, _, _, _ = ephemereal_3rules.find_pers_ephem_transients(Gcopy, persistence_thr=0.0)

        tot_graph_weight = ephemereal_3rules.compute_grossflow(Gcopy)
        graph_score = (tot_trans_weight*1.0) / (tot_graph_weight*1.0)

        print "D CLEANUP %i | Graph %s score: %f after merging nodes %s (w = %f)" % (it, filename, graph_score, bt, w)
        print "D CLEANUP %i |    merging %s | %s" % (it, u_d, v_d)
        #print "CLEANUP",it, "| trans", sorted(transients.items())
    
    it+=1
    fname = ''.join(("SM_cleanup_iter_",str(it),"_net.gexf"))
    print "<> iteration",it,"-save gexf file here:", fname
    
    nx.write_gexf(Gcopy, fname)

    dump_net_json(Gcopy, fname)


#update bad_x_transients list (may have changed)    
#Gcopy=g.copy()
transients, tot_trans_weight, bad_x_transients_full, _, edges_to_go = ephemereal_3rules.find_pers_ephem_transients(Gcopy, persistence_thr=0.0)
sorted_bad_x_transients_full = list( sorted( bad_x_transients_full.items(), key=itemgetter(1) ) )

#remove weird t2 edges
#print "EDGES_TO_GO", edges_to_go
for u,v in edges_to_go:
    print "Removing t2 edge",u,v,"because weird."
    Gcopy.remove_edge(u,v)
    

###and get rid of bad X transients

#for bt, w in sorted_bad_x_transients_full: #((focal_node, (links)), weight)
#updating transients list every time, we need to use a while contruct
it=0
while ( len(sorted_bad_x_transients_full)>0 ):
    #print "Bad X-transients pipe contains:", sorted_bad_x_transients_full

    #take out weird x transients (to exclude, they are good)
    #for (bt,w) in x_weird_transients:    
    #    sorted_bad_x_transients_full.remove((bt,w))
    #print "--- weird x transients excluded:", x_weird_transients

    bt, w = sorted_bad_x_transients_full.pop(0)
    
    focal_node = bt[0]
    bt_t2_links = bt[1]
    
        
    focal_year, focal_comm = split_node_name(focal_node)    
    
    if len(bt_t2_links)>2: #what a mess, let's put it one one side for the moment
        
        #let's treat it two at a time:
        bt_t2_links = bt_t2_links[0:2]
        print "Multi-t2 links bad X-transient; let's take two t2 links at a time:", bt_t2_links

    #we need the articles within each node/community to compute for the split. Are they anywhere to be found?
    [commsize_y0, arts_y0, G_y0] = hu.read_gexf_year(focal_year)


    print "BT: (test) focal_node, size", bt, bt[0], Gcopy.node[focal_node]['size']
    
    nodes_y0 = arts_y0[focal_comm]
    
    #nodes in stream 1 (bt_t2_links[0])
    u_0_0 = bt_t2_links[0][0]
    u_0_1 = bt_t2_links[0][1]

    u_1_0 = bt_t2_links[1][0]
    u_1_1 = bt_t2_links[1][1]
    
    #sometimes we have t2 links cross-linking 00 and 01 AND 00 and 11. Crap.
    #Let's just leve these alone (AFAF?) ##check in ephemereal library?
    if ((u_0_0 == u_1_0) or (u_0_1 == u_1_1)): #same on the left OR on the right
        print "Weird t2 link (predecessor != successor), it's not a bad transient."
        x_weird_transients.update((bt,w))
        continue

    #stream 0, year - 1 and year + 1 (in the middle, focal node): year and community
    ym1_str0_year, ym1_str0_comm = split_node_name(u_0_0)
    yp1_str0_year, yp1_str0_comm = split_node_name(u_0_1)

    #stream 1, year - 1 and year + 1 (in the middle, focal node): year and community
    ym1_str1_year, ym1_str1_comm = split_node_name(u_1_0)
    yp1_str1_year, yp1_str1_comm = split_node_name(u_1_1)

    #stream 0, extract nodes (top left and right, to help visualize)
    nodes_ym1_str0 = hu.read_gexf_year_comm(focal_year-1, ym1_str0_comm)
    nodes_yp1_str0 = hu.read_gexf_year_comm(focal_year+1, yp1_str0_comm)

    #stream 0, extract nodes (top left and right, to help visualize)
    nodes_ym1_str1 = hu.read_gexf_year_comm(focal_year-1, ym1_str1_comm)
    nodes_yp1_str1 = hu.read_gexf_year_comm(focal_year+1, yp1_str1_comm)

    #nodes from stream #0 and #1 (y-1 AND y+1)
    nodes_str0 = set(nodes_ym1_str0+nodes_yp1_str0)
    nodes_str1 = set(nodes_ym1_str1+nodes_yp1_str1)
    
    #beware: there MAY be nodes in both stream0 and stream1, let's set them aside
    focal_node_split0 = [u for u in nodes_y0 if u in nodes_str0 and u not in nodes_str1]
    focal_node_split1 = [u for u in nodes_y0 if u in nodes_str1 and u not in nodes_str0]
    
    
    focal_node_leftover = [u for u in nodes_y0 if u not in (focal_node_split0+focal_node_split1)]
    #what to do with leftovers?
    #check neighbors, assign to the half with higher links values
    for lo_u in focal_node_leftover:
        str0_acc=0
        str1_acc=0
        u_nb = G_y0.neighbors(lo_u)
        for uu_nb in u_nb:
            #count prevalent neighbors (ugly)
            if uu_nb in nodes_str0:
                str0_acc += 1
            if uu_nb in nodes_str1:
                str1_acc += 1
                
        #print "neighbor",uu_nb, str0_acc, str1_acc
        if str0_acc > str1_acc:
            focal_node_split0.append(lo_u) #move to half 0
        elif str0_acc < str1_acc:
            focal_node_split1.append(lo_u) #move to half 0
        else:
            if random()>0.5:
                focal_node_split0.append(lo_u) #move to half 0
            else:                
                focal_node_split1.append(lo_u) #move to half 0
        #we need to check that no communities are made by a single node
        if len(focal_node_split0)==1:
            #print "WARNING: very small community resulting from split. Re-run."
            #sys.exit(1)
            focal_node_split0.append(focal_node_split1.pop(0))
            
        if len(focal_node_split1)==1:
            #print "WARNING: very small community resulting from split. Re-run."
            #sys.exit(1)
            focal_node_split1.append(focal_node_split0.pop(0))
        
    size_stream0 = len(focal_node_split0)
    size_stream1 = len(focal_node_split1)
    

    #Now, split the node
    #focal_node -> focal_0 focal_1

    #rewire focal_node edges: 
    #for stream 0, link to split0, for stream 1, link to split1

    focal_node_split0_label = 's'.join(('0',focal_node)) #AFAF FIX NAMES TOO
    focal_node_split1_label = 's'.join(('1',focal_node))
    Gcopy.add_node(focal_node_split0_label, deepcopy(Gcopy.node[focal_node]))
    Gcopy.add_node(focal_node_split1_label, deepcopy(Gcopy.node[focal_node]))
    Gcopy.node[focal_node_split0_label]['label']=focal_node_split0_label
    Gcopy.node[focal_node_split1_label]['label']=focal_node_split1_label
    
    #SET CORRECT SIZE (after split)
    Gcopy.node[focal_node_split0_label]['size']=size_stream0
    Gcopy.node[focal_node_split1_label]['size']=size_stream1

    #print "old links", focal_node, u_0_0, u_0_1, u_1_0, u_1_1
    #print "focal neighs", Gcopy.neighbors(focal_node)
    ed_00 = deepcopy( Gcopy[focal_node][u_0_0] )#data for edge -> 0_0, then to link split node 0
    ed_01 = deepcopy( Gcopy[focal_node][u_0_1] ) # 0_1
    ed_10 = deepcopy( Gcopy[focal_node][u_1_0] ) # 1_0
    ed_11 = deepcopy( Gcopy[focal_node][u_1_1] ) # 1_1
    #fix ids (trash them? no, be gentle and relabel)
    ed_00_id = 's'.join(('0', ed_00['id']))
    ed_01_id = 's'.join(('1', ed_00['id']))
    ed_10_id = 's'.join(('2', ed_00['id']))
    ed_11_id = 's'.join(('3', ed_00['id']))

    ed_00['id'] = ed_00_id
    ed_01['id'] = ed_01_id
    ed_10['id'] = ed_10_id
    ed_11['id'] = ed_11_id


    #AFAFAF FIX Jaccards!
    #relink new top half
    #if focal_node=='3_1995':
    #    print "3_1995***", Gcopy.neighbors(focal_node)
    Gcopy.add_edge(focal_node_split0_label, u_0_0, ed_00)
    Gcopy.add_edge(focal_node_split0_label, u_0_1, ed_01)
    Gcopy.add_edge(focal_node_split1_label, u_1_0, ed_10)
    Gcopy.add_edge(focal_node_split1_label, u_1_1, ed_11)

    #if focal_node=='3_1995':
    #    print "3_1995_s0***", Gcopy.neighbors(focal_node_split0_label)
    #    print "3_1995_s1***", Gcopy.neighbors(focal_node_split1_label)
    
    Gcopy.remove_node(focal_node)


    #update history master 
    #[0] - {comm: size}
    #[1] - {comm: arts list}
    #[2] - {G}
    #-> nodes split
    
    #[1] fix resulting merged community
    new_comm_label0 = change_comm_name_split(focal_comm, 0)
    new_comm_label1 = change_comm_name_split(focal_comm, 1)
    print "UUU", focal_year, focal_comm, new_comm_label0, new_comm_label1, nodes_y0, size_stream0, size_stream1, focal_node_split0, focal_node_split1,\
        focal_node_leftover
        
    if focal_comm=='0':
        print "cn0:", new_comm_label0, new_comm_label1

    #[0]
    history_master[focal_year][0][new_comm_label0] = size_stream0
    history_master[focal_year][0][new_comm_label1] = size_stream1
    #print "1deleting", focal_year, focal_comm
    #print "1...", history_master[focal_year][0]
    del history_master[focal_year][0][focal_comm]
    
    #[1]
    history_master[focal_year][1][new_comm_label0] = focal_node_split0
    history_master[focal_year][1][new_comm_label1] = focal_node_split1

    
    if focal_comm not in history_master[focal_year][1]:
        del history_master[focal_year][1][change_comm_name_merge(focal_comm)]
        #nodes_from_mn = history_master[mn_year][1][change_comm_name_merge(mn_comm)]
    else:
        del history_master[focal_year][1][focal_comm]
        #nodes_from_mn = history_master[mn_year][1][mn_comm]

    #########del history_master[focal_year][1][focal_comm]
    #[2]    
    #history_master[focal_year][2] = Gcopy
    #print "XXX", mn_year, "XXX", Gcopy.nodes()

    
    #update
    transients, tot_trans_weight, bad_x_transients_full, _, edges_to_go = ephemereal_3rules.find_pers_ephem_transients(Gcopy, persistence_thr=0.0)
    #remove weird t2 edges
    #print "EDGES_TO_GO", edges_to_go
    for u,v in edges_to_go:
        print "Removing t2 edge",u,v,"because weird."
        Gcopy.remove_edge(u,v)



    sorted_bad_x_transients_full = list( sorted( bad_x_transients_full.items(), key=itemgetter(1) ) )

    tot_graph_weight = ephemereal_3rules.compute_grossflow(Gcopy)
    graph_score = (tot_trans_weight*1.0) / (tot_graph_weight*1.0)

    it+=1
    print "X CLEANUP %i | Graph %s score: %f after splitting node %s (w = %f)" % (it, filename, graph_score, bt[0], w)
    print "X CLEANUP    |    s0", focal_node_split0
    print "X CLEANUP    |    s1", focal_node_split1
    #print "X CLEANUP %i |    merging %s | %s" % (it, u_d, v_d)

    fname = ''.join(("MS_cleanup_iter_",str(it),"_net.gexf"))
    print "X iteration",it,"-save gexf file here:", fname
    
    
    ################################################################MMM
    """
    if focal_node == '0_1994':
        print "<<<>>>", size_stream0, focal_node_split0_label
        print "<<<>>>", focal_node_split0
        print "<<<>>>", size_stream1, focal_node_split1_label
        print "<<<>>>", focal_node_split1
    """ 

    nx.write_gexf(Gcopy, fname)

    dump_net_json(Gcopy, fname)

    #with open("history_master.json", "wb") as hmoutfile:
    #    json.dump(history_master, hmoutfile)

#sys.exit(0)

#COLOR TRANSIENTS FOR FINAL 
#ignore full transients structure (args 2,3)
transients, tot_trans_weight, _, _, _= ephemereal_3rules.find_pers_ephem_transients(Gcopy, persistence_thr=0.0) 
tot_graph_weight = ephemereal_3rules.compute_grossflow(Gcopy)

graph_score = (tot_trans_weight*1.0) / (tot_graph_weight*1.0)

print "Final graph score: %f" % (graph_score)
#print "Untreated transients:", untreated_transients
print "Weird X transients:", x_weird_transients

#now, refresh idcards with modified communities
#def idcard_year_splitmerge(year, partitions_size, list_nodes, G):
for y in years:
    print "Making updated IDCard for year", y
    [partitions_size, list_nodes, G] = history_master[y]
    #idcm.idcard_year_splitmerge(y, partitions_size, list_nodes, G)

for u,d in Gcopy.nodes(data=True):
    #pos
    u_x=d['Xpos']
    u_y=d['y-pos-fix']
    
    Gcopy.node[u]['viz']={}
    Gcopy.node[u]['viz']['position']={}

    Gcopy.node[u]['viz']['position']['x'] = u_x*100
    Gcopy.node[u]['viz']['position']['y'] = u_y*10
    Gcopy.node[u]['viz']['position']['z'] = 0.0

    #initialize spurious transients colors
    Gcopy.node[u]['trans_col'] = 0

for tr, w in transients.items():
    if w>0:        
        if isinstance(tr, basestring): #X-transients, one node, detect if is the only element
            Gcopy.node[tr]['trans_col'] = 1
            #print "Coloring node",tr,"as good (x-trans)"
        else: #<>-transients, more than one node
            for u in tr:
                Gcopy.node[u]['trans_col'] = 1
                #print "Coloring node",u,"as good (<>-trans)"
    elif w<0: #there are no zeros, but better safe than sorry
        if isinstance(tr, basestring): #X-transients, one node
            Gcopy.node[tr]['trans_col'] = 2
            #print "Coloring node",tr,"as bad (x-trans)"
        else:
            for u in tr:
                Gcopy.node[u]['trans_col'] = 2
                #print "Coloring node",u,"as bad (<>-trans)"

####compute scores (skifo, modularity) for POST GRAPHS
for y in years:
    [_, _part, _G] = history_master[y] 
    _part_as_nodes_comms = part2nodes(_part)
    #some nodes are in communities filtered out because too small, take them out of G when computing modularity
    lost_nodes = [u for u in _G.nodes() if u not in _part_as_nodes_comms.keys()]
    _G_clean = _G.copy()
    for u in lost_nodes:
        _G_clean.remove_node(u)
    
    #print "ln", len(lost_nodes), lost_nodes
    _mod = community.modularity(_part_as_nodes_comms, _G_clean)
    #_mod = community.modularity(_part_as_nodes_comms, _G)
    print "POST-SCORES", y, _mod
                
nx.write_gexf(Gcopy, 'net_final.gexf')

json_filename = 'net_final'
dump_net_json(Gcopy, json_filename)

print "DUMPING FINAL NODES by PARTITIONS"

for y in years:
    [_, _part, _] = history_master[y] 

    json_filename = 'part_'+str(y)+'.json'
    with open(json_filename, 'w') as outfile:
        json.dump(_part, outfile)


"""
#compute filtered nodes weight
m_nodes=['0_2003', '26_2003', '15_2003', '64_2003', '0_1994', '12_1994',
'10_1997', '21_1997',
'0_2000', '1_2000',
'11_2005', '40_2005',
'14_2005', '57_2005',
'85_2003', '86_2003',
'52_1997', '43_1997',
'9_1989', '8_1989',
'40_1986', '46_1986',
'19_1981', '4_1981']

s_nodes=['0_2003',
'0_1994', 
'1_1999' ,
'7_2001' ,
'17_2005',
'3_1992' ,
'3_2008',
'10_2006',
'5_2001',
'12_2005',
'18_1999',
'20_2002']

ms_f_nodes=set(m_nodes+s_nodes)

w_all_nodes = sum([ k['size'] for u,k in g.nodes(data=True) ])
w_m_nodes = sum([ k['size'] for u,k in g.nodes(data=True) if u in m_nodes ])
w_s_nodes = sum([ k['size'] for u,k in g.nodes(data=True) if u in s_nodes ])
w_ms_f_nodes = sum([ k['size'] for u,k in g.nodes(data=True) if u in ms_f_nodes ])



print "Total weight:", w_all_nodes
print "Merged nodes weight:", w_m_nodes, (w_m_nodes*1.0/w_all_nodes)*100.0
print "Split nodes weight:", w_s_nodes, (w_s_nodes*1.0/w_all_nodes)*100.0

print "Split+merged nodes weight:", w_ms_f_nodes, (w_ms_f_nodes*1.0/w_all_nodes)*100.0
"""