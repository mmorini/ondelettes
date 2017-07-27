#read json community-artifcles files, measures jaccard index
#parts (C) SG

import os
#import sys
#import glob
#import numpy
#import Utils
#import time
#import json
#import math
import networkx as nx
#import matplotlib.pyplot as plt
#import operator
import histUtils as hu
#AFAF : throw away weighted graph??

#import csv

datadir='/home/matteo/WORK/PhD/projets/ONDELETTES/FINAL/1.2.3'

firstyear=1970
lastyear=2010

comm_size_threshold = 5

def main():
    
    #G=nx.DiGraph() #unweighted (one link only)
    Gw=nx.DiGraph() #weighted (one predecessor, one successor)

    
    for y in range(firstyear,lastyear-2): #range goes up to the penultimate, 2008-1= up to 2007(+8+9)
        y1=y
        y2=y+1
        y3=y+2
    
        #both unweighted and weighted edges (all) are returned
        print "COMPUTING JACCARD on years",y1,y2
        [nodes1, nodes2, nodes3, edges, keyword_dict1, auth_dict1,\
            keyword_dict2, auth_dict2, keyword_dict3, auth_dict3] = jaccard(y1,y2, thr = comm_size_threshold) 
        """
        print "****************"
        print nodes1
        print nodes2
        print edges
        print "****************"
        """
        Gw.add_nodes_from(nodes1)
        Gw.add_nodes_from(nodes2)
        Gw.add_nodes_from(nodes3)
        Gw.add_edges_from(edges)
            
        _nodes=Gw.nodes()
        for n in _nodes: #same loop for Gw, profiting

            if (n in nodes1):
                #G.node[n]['year'] = y1
                #G.node[n]['size'] = nodes1[n]
                #G.node[n]['Xpos'] = y1-firstyear

                Gw.node[n]['year'] = y1
                Gw.node[n]['size'] = nodes1[n]
                Gw.node[n]['Xpos'] = y1-firstyear

                c1_noyear=int(n.split("_")[0])                
                if c1_noyear in keyword_dict1:
                #    G.node[n]['keyword'] = keyword_dict1[c1_noyear]
                    Gw.node[n]['keyword'] = keyword_dict1[c1_noyear]
                else:
                #    G.node[n]['keyword'] = "N/A"
                    Gw.node[n]['keyword'] = "N/A"
                    
                #G.node[n]['mainauth'] = auth_dict1[c1_noyear]
                Gw.node[n]['mainauth'] = auth_dict1[c1_noyear]
                #pseudogeo layout
                #G.node[n]['x'] = ((y1-1900.0)*10.0)
                Gw.node[n]['x'] = ((y1-1900.0)*10.0)
                if c1_noyear>999:
                    c1_ny_corr=(c1_noyear/1000.0)+100
                else:
                    c1_ny_corr=c1_noyear
                #G.node[n]['y'] = c1_ny_corr*1.0
                Gw.node[n]['y'] = c1_ny_corr*1.0
                
                
                #G.node[n]['viz']=dict()
                ##G.node[n]['viz']['position']=dict()
                #G.node[n]['viz']['position']['x'] = G.node[n]['x']
                #G.node[n]['viz']['position']['y'] = G.node[n]['y']
                #G.node[n]['viz']['position']['z'] = 0.0
                
                Gw.node[n]['viz']=dict()
                Gw.node[n]['viz']['position']=dict()
                Gw.node[n]['viz']['position']['x'] = Gw.node[n]['x']
                Gw.node[n]['viz']['position']['y'] = Gw.node[n]['y']
                Gw.node[n]['viz']['position']['z'] = 0.0
                

            elif (n in nodes2):
                #G.node[n]['year'] = y2
                #G.node[n]['size'] = nodes2[n]
                #G.node[n]['Xpos'] = y2-firstyear


                Gw.node[n]['year'] = y2
                Gw.node[n]['size'] = nodes2[n]
                Gw.node[n]['Xpos'] = y2-firstyear

                c2_noyear=int(n.split("_")[0])                
                if c2_noyear in keyword_dict2:
                #    G.node[n]['keyword'] = keyword_dict2[c2_noyear]
                    Gw.node[n]['keyword'] = keyword_dict2[c2_noyear]
                else:
                #    G.node[n]['keyword'] = "N/A"
                    Gw.node[n]['keyword'] = "N/A"

                #print "auth_dict2", auth_dict2
                #G.node[n]['mainauth'] = auth_dict2[c2_noyear]
                Gw.node[n]['mainauth'] = auth_dict2[c2_noyear]
                #pseudogeo layout
                #G.node[n]['x'] = ((y2-1900.0)*10.0)
                Gw.node[n]['x'] = ((y2-1900.0)*10.0)
                if c2_noyear>999:
                    c2_ny_corr=(c2_noyear/1000.0)+100
                else:
                    c2_ny_corr=c2_noyear
                #G.node[n]['y'] = c2_ny_corr*1.0            
                Gw.node[n]['y'] = c2_ny_corr*1.0            
    
    
                #G.node[n]['viz']=dict()
                #G.node[n]['viz']['position']=dict()
                #G.node[n]['viz']['position']['x'] = G.node[n]['x']
                #G.node[n]['viz']['position']['y'] = G.node[n]['y']
                #G.node[n]['viz']['position']['z'] = 0.0
                
                Gw.node[n]['viz']=dict()
                Gw.node[n]['viz']['position']=dict()
                Gw.node[n]['viz']['position']['x'] = Gw.node[n]['x']
                Gw.node[n]['viz']['position']['y'] = Gw.node[n]['y']
                Gw.node[n]['viz']['position']['z'] = 0.0
                
            elif (n in nodes3):
                #G.node[n]['year'] = y3
                #G.node[n]['size'] = nodes3[n]
                #G.node[n]['Xpos'] = y3-firstyear

                Gw.node[n]['year'] = y3
                Gw.node[n]['size'] = nodes3[n]
                Gw.node[n]['Xpos'] = y3-firstyear

                c3_noyear=int(n.split("_")[0])                
                if c3_noyear in keyword_dict3:
                #    G.node[n]['keyword'] = keyword_dict3[c3_noyear]
                    Gw.node[n]['keyword'] = keyword_dict3[c3_noyear]
                else:
                #    G.node[n]['keyword'] = "N/A"
                    Gw.node[n]['keyword'] = "N/A"

                #print "auth_dict2", auth_dict2
                #G.node[n]['mainauth'] = auth_dict3[c3_noyear]
                Gw.node[n]['mainauth'] = auth_dict3[c3_noyear]
                #pseudogeo layout
                #G.node[n]['x'] = ((y3-1900.0)*10.0)
                Gw.node[n]['x'] = ((y3-1900.0)*10.0)
                if c3_noyear>999:
                    c3_ny_corr=(c3_noyear/1000.0)+100
                else:
                    c3_ny_corr=c3_noyear
                #G.node[n]['y'] = c3_ny_corr*1.0            
                Gw.node[n]['y'] = c3_ny_corr*1.0            
    
                #G.node[n]['viz']=dict()
                #G.node[n]['viz']['position']=dict()
                #G.node[n]['viz']['position']['x'] = G.node[n]['x']
                #G.node[n]['viz']['position']['y'] = G.node[n]['y']
                #G.node[n]['viz']['position']['z'] = 0.0

                Gw.node[n]['viz']=dict()
                Gw.node[n]['viz']['position']=dict()
                Gw.node[n]['viz']['position']['x'] = Gw.node[n]['x']
                Gw.node[n]['viz']['position']['y'] = Gw.node[n]['y']
                Gw.node[n]['viz']['position']['z'] = 0.0
                
    #nx.write_gexf(G, "./g_big.gexf", prettyprint=True)
    nx.write_gexf(Gw, "./g_big_w.gexf", prettyprint=True)

    
def jaccard(year1,year2, thr = 5): #and year3 (y2+1); default threshold = 5 (min comm size)
    
    y1=year1

    y2=year2
    y2last=y2+3

    y3=year2+1 #t+2
    y3last=y3+3 #t+2
    
    d_y1=''.join(['data',str(y1)])
    d_y2=''.join(['data',str(y2)])
    d_y3=''.join(['data',str(y3)])
    #d_y2last=''.join(['data',str(y2last)])
    
    dir1 = os.path.join(datadir,str(d_y1))
    dir2 = os.path.join(datadir,str(d_y2))
    dir3 = os.path.join(datadir,str(d_y3))
    
    #filename='partitions_thr2.txt'
    filename_partitions='partitions.gexf'
    filename_nodes='nodes.gexf'
        
    gexf_filename1 = os.path.join(dir1, filename_nodes)
    gexf_filename2 = os.path.join(dir2, filename_nodes)
    gexf_filename3 = os.path.join(dir3, filename_nodes)
    part_filename1 = os.path.join(dir1, filename_partitions)
    part_filename2 = os.path.join(dir2, filename_partitions)
    part_filename3 = os.path.join(dir3, filename_partitions)

    #print "LOADING FILES", gexf_filename1, gexf_filename2, part_filename1, part_filename2
    
    #read json files
    
    #rather, load GEXF files!
    [commsize1, arts1, G1] = hu.read_gexf_net(gexf_filename1, thr)
    [commsize2, arts2, G2] = hu.read_gexf_net(gexf_filename2, thr)
    [commsize3, arts3, G3] = hu.read_gexf_net(gexf_filename3, thr) #t+2 too

    #prepare author, keyword dictionaries    
    G1_part = nx.read_gexf(part_filename1)
    G2_part = nx.read_gexf(part_filename2)
    G3_part = nx.read_gexf(part_filename3)


    #print "COMMSIZE1: ", y1, len(commsize1)

    #update nodes names: comms and arts
    keyword1_dict={}
    author1_dict={}
    keyword2_dict={}
    author2_dict={}
    keyword3_dict={}
    author3_dict={}
    
    k1=commsize1.keys()   
    for c in k1:
        newkey=str(c)+"_"+str(y1)
        #print commsize1, c
        commsize1[newkey]=commsize1.pop(c)      
        keyword1_dict[int(c)]=G1_part.node[c]['comm_label3']
        author1_dict[int(c)]=G1_part.node[c]['main_author']
        
    k2=commsize2.keys()       
    for c in k2:
        newkey=str(c)+"_"+str(y2)
        commsize2[newkey]=commsize2.pop(c)
        keyword2_dict[int(c)]=G2_part.node[c]['comm_label3']
        author2_dict[int(c)]=G2_part.node[c]['main_author']

    k3=commsize3.keys()       
    for c in k3:
        newkey=str(c)+"_"+str(y3)
        commsize3[newkey]=commsize3.pop(c)
        keyword3_dict[int(c)]=G3_part.node[c]['comm_label3']
        author3_dict[int(c)]=G3_part.node[c]['main_author']
        
    #relabel nodes, communities get year appended as in xx_year
    k1a=arts1.keys()
    for c in k1a:
        newkey=str(c)+"_"+str(y1)
        #print "commsize1, c", commsize1, c
        arts1[newkey]=arts1.pop(c)
        
    k2a=arts2.keys()
    for c in k2a:
        #print "c", k2a, c
        newkey=str(c)+"_"+str(y2)
        #print commsize1, c
        arts2[newkey]=arts2.pop(c)

    k3a=arts3.keys()
    for c in k3a:
        #print "c", k2a, c
        newkey=str(c)+"_"+str(y3)
        #print commsize1, c
        arts3[newkey]=arts3.pop(c)
    
    ## COMPUTE SIMILARITIES - ARTICLES
    print "Jaccard time"  

    jac = dict()
    jac_t2 = dict()

    #let's clean up this mess, please AFAF
    #t+1
    for c1 in commsize1: #1-2
        
        subset_arts1 = [art for art in arts1[c1] if G1.node[art]['year'] > y1] #remove first year

        for c2 in commsize2:

            subset_arts2 = [artb for artb in arts2[c2] if G2.node[artb]['year'] < y2last] #remove last year
            
            shared_art=[art for art in subset_arts1 if art in subset_arts2]
            onlyin1_art=[art for art in subset_arts1 if art not in subset_arts2]
            onlyin2_art=[art for art in subset_arts2 if art not in subset_arts1]
            
            #.. jaccard index: |intersection| / |union|
            #c1,c2 must then be offset by one to match BC pdfs
            jac[(c1,c2)] = len(shared_art)*1.0/(len(onlyin1_art)+len(onlyin2_art)+len(shared_art))

    for c2 in commsize2: #2-3
        #beware of slices overlaps, here 2 years only
        subset_arts2 = [art for art in arts2[c2] if G2.node[art]['year'] > y2] #remove first year

        for c3 in commsize3:

            subset_arts3 = [artc for artc in arts3[c3] if G3.node[artc]['year'] < y3last] #remove last year
            
            shared_art=[art for art in subset_arts2 if art in subset_arts3]
            onlyin2_art=[art for art in subset_arts2 if art not in subset_arts3]
            onlyin3_art=[art for art in subset_arts3 if art not in subset_arts2]
            #.. jaccard index: |intersection| / |union|
            #c1,c2 must then be offset by one to match BC pdfs
            jac[(c2,c3)] = len(shared_art)*1.0/(len(onlyin2_art)+len(onlyin3_art)+len(shared_art))

    #t+2
    for c1 in commsize1: #1-3
        #beware of slices overlaps, here 2 years only
        subset_arts1 = [art for art in arts1[c1] if G1.node[art]['year'] > y1+1] #remove first two years

        for c3 in commsize3:

            subset_arts3 = [artc for artc in arts3[c3] if G3.node[artc]['year'] < y3last-1] #remove last two years
            
            shared_art=[art for art in subset_arts1 if art in subset_arts3]
            onlyin1_art=[art for art in subset_arts1 if art not in subset_arts3]
            onlyin3_art=[art for art in subset_arts3 if art not in subset_arts1]
            #.. jaccard index: |intersection| / |union|
            #c1,c2 must then be offset by one to match BC pdfs
            print "ppp1",c1,c3,shared_art, onlyin1_art, onlyin3_art, shared_art
            print "ppp2",c1,c3,len(shared_art), len(onlyin1_art), len(onlyin3_art), len(shared_art)
            _num = len(shared_art)*1.0
            _den = (len(onlyin1_art)+len(onlyin3_art)+len(shared_art))
            if _den != 0:
                jac_t2[(c1,c3)] = _num / _den
            else:
                jac_t2[(c1,c3)] = 0.0

    print "JAC",jac

    #precedessors and successors
    
    #successors: loop each t1 community around t2 communities and pick the highst JAC
    successors={}       #originating node : (successor, weight in jac)
    successors_all={}    
    successors_t2={} #also unique pred/succ for t+2 links
    
    for c1 in commsize1.keys(): #1->2 1->3
        succ=-1
        maxjac=0
        succ_t2=-1
        maxjac_t2=0
        successors_all[c1]=list()
        
        for c2 in commsize2.keys(): #t1
            
            if jac[(c1,c2)]>maxjac:
                maxjac=jac[(c1,c2)]
                succ=c2

            ##also 2nd, 3rd... ranked successors
            if jac[(c1,c2)]>0:
                successors_all[c1].append((c2,jac[(c1,c2)])) 

            if succ!=-1:
                successors[c1]=(succ, maxjac) #usual offset stuff
        
        for c3 in commsize3.keys(): #t2
            if jac_t2[(c1,c3)]>maxjac_t2: #!!!!!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                maxjac_t2=jac_t2[(c1,c3)]
                succ_t2=c3

            if succ_t2!=-1:
                successors_t2[c1]=(succ_t2, maxjac_t2) #usual offset stuff

    ###precedessors: loop each t2 community around t1 communities and pick the highst JAC
    predecessors={}
    predecessors_all={}    
    predecessors_t2={}

    for c3 in commsize3.keys(): #3->2 3->1   t+1 - second year looped first
        pred=-1
        maxjac=0
        pred_t2=-1
        maxjac_t2=0
        predecessors_all[c3]=list()

        for c2 in commsize2.keys():
            #print c2,c3
            if jac[(c2,c3)]>maxjac:                
                maxjac=jac[(c2,c3)]
                pred=c2
                
            ##also 2nd, 3rd... ranked successors
            if jac[(c2,c3)]>0:
                predecessors_all[c3].append((c3,jac[(c2,c3)])) 

            if pred!=-1:
                predecessors[c3]=(pred, maxjac)


        for c1 in commsize1.keys():
            #print jac[(c1,c2)]
            
            if jac_t2[(c1,c3)]>maxjac_t2: 
                maxjac_t2=jac_t2[(c1,c3)]
                pred_t2=c1
                
            if pred_t2!=-1:
                predecessors_t2[c3]=(pred_t2, maxjac_t2) #usual offset stuff

    ###################### edges
    edges=list()
    #t+1's
    for p, d in predecessors.items():
        #print p,d
        e0 = p
        e1 = d[0]
        w = d[1]        
        edges.append( (e0, e1, {'weight': w, 'scope': 't1'} ) )

    for s, d in successors.items():
        #print s,d
        e0 = s
        e1 = d[0]
        w = d[1]        
        edges.append( (e0, e1, {'weight': w, 'scope': 't1'} ) )

    #t+2's
    for p, d in predecessors_t2.items():
        #print p,d
        e0 = p
        e1 = d[0]
        w = d[1]
        edges.append( (e0, e1, {'weight': w, 'scope': 't2'} ) )

    for s, d in successors_t2.items():
        #print s,d
        e0 = s
        e1 = d[0]
        w = d[1]        
        edges.append( (e0, e1, {'weight': w, 'scope': 't2'} ) )

    #also return keywords and authors as dictionaries

    return [commsize1, commsize2, commsize3, edges,         keyword1_dict, author1_dict, keyword2_dict, author2_dict, keyword3_dict, author3_dict]
    


if __name__ == "__main__":
    main()
    
