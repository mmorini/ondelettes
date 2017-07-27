import os
from networkx import nx

def dump_json_nodescomms():
    nodesc
    

def read_gexf_net(file, thr = 5): #minimum community size = 5
    #lev ignored (only 'atomic' level present)
    if not os.path.isfile(file):
        print '....file %s does not exist' % file
        exit(-1)
    else:
        print "read_gexf loading file", file 

    G = nx.read_gexf(file)
    
    partitions = set([d['part'] for u,d in G.nodes(data=True)])
    
    list_nodes = dict();
    comm_size = dict();

    for com in set(partitions) :

        list_nodes[com] = [node for node, data in G.nodes(data = True) if data['part'] == com]
        size = len(list_nodes[com])
    
        #only communities at least n sized
        #print "size", size
        if size > thr:
            comm_size[com] = size

    return [comm_size, list_nodes, G]

def read_gexf_year(year, thr = 5):

    datadir='/home/matteo/WORK/PhD/projets/ONDELETTES/FINAL/1.2.3'
    filename_nodes='nodes.gexf'
    d_y0       = ''.join(['data',str(year)])
    _dir_y0     = os.path.join(datadir,str(d_y0))
    _filename_y0     = os.path.join(_dir_y0, filename_nodes)

    [commsize, arts, G] = read_gexf_net(_filename_y0)

    return [commsize, arts, G]

def read_gexf_year_comm(year, comm, thr = 5):
    
    [commsize, arts, G] = read_gexf_year(year)
    
    #print "...", arts
    nodes = arts[comm]
    
    return nodes
