# -*- coding: utf-8 -*-
"""
Created on Fri May 19 10:32:43 2017

@author: matteo
"""

#IDCards maker
import os
import sys
import glob
import numpy
import argparse
import time
import json
import math
import Utils, BCUtils
import networkx as nx
#import community_atoms
#from subprocess import call, Popen, PIPE
import histUtils as hu


##########################################################################################
#prepare data
#y1=min([d['year'] for u,d in G.nodes(data=True)])
#y2=max([d['year'] for u,d in G.nodes(data=True)])

#year = 1990


############################################à


def main():
    print "Preparing .tex IDCards files."
    idcards_allyears()
    print "Done."

    
#first, original graph:
#years 1971-2009

#input: year

def idcards_allyears():
    for year in range(1971,2009+1):
        print "...working on year", year
        idcard_year(year)
        print "...exiting id_year function."
    return


def idcard_year(year):

    #partitions_size filters only comms bigger than thr
    #partitions_size used to pick big enough partitions, nodes contains all nodes, G
    [partitions_size, list_nodes, G] = hu.read_gexf_year(year) #list_nodes = {comm: [list nodes]}
    idcard_year_splitmerge(year, partitions_size, list_nodes, G)
        
def idcard_year_splitmerge(year, partitions_size, list_nodes, G):

    print "idcard_year_splitmerge | year", year
    verbose = True
    in_dir = '/home/matteo/WORK/PhD/projets/ONDELETTES/FINAL/1.2.3/dat_files_clean'
    out_dir_root = '/home/matteo/WORK/PhD/projets/ONDELETTES/FINAL/6'
    
    thr = 5 #AFAFfix
    bcthr = thr
    
    txt = 'level_atom' 
    textxt = 'level atom' 
    
    #Gfilename = 'net_final_layedout_greens_only.gexf'
    #Gfilename = 'g_big_w.gexf_viz.gexf'
    #G = nx.read_gexf(Gfilename)
    
    #duplicated later, fix, AFAF
    src1  = os.path.join(in_dir, "articles.dat") 
    pl = Utils.Article()
    pl.read_file(src1)
    nb_art = len(pl.articles) # store the number of articles within database

    out_dir = os.path.join(out_dir_root,str(year))
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fooname="BCcomm_ID_Cards_thr%d_%s.tex" % (bcthr,txt)
    filename = os.path.join(out_dir, fooname)
    f_out = open(filename,"w")
    f_out.write("\documentclass[a4paper,11pt]{report}\n\usepackage[english]{babel}\n\usepackage[latin1]{inputenc}\n\usepackage{amsfonts,amssymb,amsmath}\n\usepackage{pdflscape}\n\usepackage{color}\n\usepackage{graphicx}\n\n\\addtolength{\evensidemargin}{-60pt}\n\\addtolength{\oddsidemargin}{-60pt}\n\\addtolength{\\textheight}{80pt}\n\n\\title{{\\bf Communities ID Cards}}\n\date{\\begin{flushleft}This document gathers the ``ID Cards'' of the BC communities found within the studied database.\\\\\n The BC network was built by linking pairs of articles sharing at least %d references - %d out %d articles are in the network. The %d communities presented here correspond to the ones found in the %s level grouping at least %d articles. They correspond to a total of %d articles. \\\\\n These ID cards displays the most frequent keywords, subject categories, journals of publication, institution, countries, authors, references and reference journals of the articles of each community. The significance of an item $\sigma = \sqrt{N} (f - p) / \sqrt{p(1-p)}$ [where $N$ is the number of articles within the community and $f$ and $p$ are the proportion of articles respectively within the community and within the database displaying that item ] is also given (for example $\sigma > 5$ is really highly significant). Each community is labelled by the most frequent keyword with positive $\sigma$.\\\\\n\\vspace{1cm}\n\copyright Sebastian Grauwin - BiblioTools3.0 (September 2014) \end{flushleft}}\n\n\\begin{document}\n\\begin{landscape}\n\maketitle\n" % (bcthr, len(G.nodes()), nb_art, 999, textxt, thr, 999))
    
    f_out.write("\n\n\clearpage\\begin{figure}[h!]\n\\begin{center}\n%%\includegraphics[width=1.3\\textwidth]{xxx}\n\caption{{\\bf Top communities network.} The sizes of the nodes correspond to the number of articles in each community. Their colors correspond to the top community they belong to and the labels correspond to those reported in the tables below. The edges reflects shared references between communities. The thicker an edge between two communities, the more references they share. The arrows reflect a `who feed who' relationship: an arrow pointing from A to B indicates that the references shared by A and B are more similar to the references from A.}\n\end{center}\n\end{figure}\n\n")
    
    top_com = {k:v for k,v in list_nodes.items() if k in partitions_size.keys()}
    
    top_part = top_com
    
    #partition (idiotic): node: community
    partition = {}
    for c, us in list_nodes.items():
        for u in us:
            partition[int(u)] = c
    
    
    # inner weight of top com
    list_top = dict()
    topcomm_innerw = dict(); topcomm_size = dict()
    

    #for topcom in set(top_part.values()) : #partition numbers
    for topcom in top_part.keys() : #partition numbers
      #list_top[topcom] = [nodes for nodes in top_part.keys() if top_part[nodes] == topcom] #nodes in each partition
      list_top[topcom] = list_nodes[topcom]
      size = len(list_top[topcom])
      #print topcom, size
      W = 0;
      for id1 in list_top[topcom]:
        for id2 in list_top[topcom]:
          if id2 > id1 and id2 in G.edge[id1]: 
            W += G.edge[id1][id2]['weight']
      print "pppppp", size, topcom, list_top[topcom]
      W *= 2.0 / (size * (size -1))
      topcomm_innerw[topcom]= 1.0 / W
      topcomm_size[topcom] = [size,1]
    
    #
    
    #... quantitative stuff about partition
    comm_innerw = dict(); comm_size = dict(); 
    for com in top_com.keys():
      size = len(list_nodes[com])
      W = 0;
      for id1 in list_nodes[com]:
        for id2 in list_nodes[com]:
          if id2 > id1 and id2 in G.edge[id1]: 
            W += G.edge[id1][id2]['weight']
      W *= 2.0 / (size * (size -1))
      comm_innerw[com]= 1.0 / W
      #comm_size[com] = [size,topcomm_innerw[top_com[com]],topcomm_size[top_com[com]]]
      comm_size[com] = [size,topcomm_innerw[com],topcomm_size[com]]
    #...order communities by innerw (or size, depending on previously defined parameter "which") of their top and then by size 
    Lcomm_size = comm_size.items()
    
    #if which =='s': Lcomm_size.sort(cmpval)
    #if which =='d': Lcomm_size.sort(cmpvalb)
    Lcomm_size.sort(cmpval) #which = 's', sort by size
    
    #.. frequency / significance of keywords, etc...
    comm_label = dict();
    author_mostfreq='xxx'            
    #and there are:
    #stuffK - keywords, 
    #stuffTK - titlewords
    #stuffS - subjects
    #stuffJ - journals
    #stuffA - authors ..., stuffI, stuffC, stuffR, stuffRJ) = BCUtils.comm_tables(in_dir,partition,thr,verbose)
    (stuffK, stuffTK, stuffS, stuffJ, stuffA, stuffI, stuffC, stuffR, stuffRJ) = BCUtils.comm_tables(in_dir,partition,thr,verbose)
    
    src1  = os.path.join(in_dir, "articles.dat") 
    pl = Utils.Article()
    pl.read_file(src1)
    nb_art = len(pl.articles) # store the number of articles within database
    
    #.. output tables
    for elm in Lcomm_size:
      if elm[1][0] > thr:
        com = elm[0]
    
        if com in stuffTK:
            author_mostfreq=stuffA[com][0]
    
        #fill perhaps missingvalues (early years)
        comm_label[com]=['xxx','xxx','xxx','xxx']
    
    
        #K
        if com in stuffK:
        #if (com in stuffK) or (com in stuffTK):
    
          if len(stuffK[com]) > 0 : 
            # the different labels are based on a combination of most freq / sign keyword
            comm_label[com] = extract_labels(stuffK[com]) 
          else: 
            # if no keyword, the different labels are based on a combination of most freq / sign title word  
            comm_label[com] = extract_labels(stuffTK[com]) 
          #    
    
          #### community : subject,most freq author, size
          #ids[com]=(comm_label[com][0], comm_label[com][3], author_mostfreq[0], comm_size[com][0])
          #print "YYY", (comm_label[com][3], com, comm_size[com][0], comm_innerw[com] ) 
          f_out.write("\clearpage\n\n\\begin{table}[!ht]\n\caption{The community ``%s'' (id %d) contains $N = %d$ articles. Its average internal link weight is $<\omega_{in}> \simeq 1/%d$ }\n\\textcolor{white}{aa}\\\\\n{\scriptsize\\begin{tabular}{|l r r|}\n\hline\nKeyword & f(\\%%) & $\sigma$\\\\\n\hline\n" % (comm_label[com][3], int(com), int(comm_size[com][0]), comm_innerw[com] ) )
          for i in range(len(stuffK[com])):
            if len(stuffK[com][i][0]) < 30:
              f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffK[com][i][0], stuffK[com][i][1], stuffK[com][i][2]) )
            else:
              aux = stuffK[com][i][0].rfind(' ')
              while aux > 30: 
                 aux = stuffK[com][i][0][0:aux].rfind(' ')
              f_out.write("%s &  & \\\\\n" % ( stuffK[com][i][0][0:aux] ) )
              f_out.write("$\quad$%s & %1.2f & %1.2f\\\\\n" % ( stuffK[com][i][0][aux:], stuffK[com][i][1], stuffK[com][i][2]) )
          for i in range(max(0,20-len(stuffK[com]))):
            f_out.write(" &  & \\\\\n")
        else:
          comm_label[com] = ['xxx','xxx','xxx','xxx']
          #print "XXX", (com, comm_size[com][0], comm_innerw[com] )
          f_out.write("\clearpage\n\n\\begin{table}[!ht]\n\caption{The community ``?'' (id %d) contains $N = %d$ articles. Its average internal link weight is $<\omega_{in}> \simeq 1/%d$ }\n\\textcolor{white}{aa}\\\\\n{\scriptsize\\begin{tabular}{|l r r|}\n\hline\nKeyword & f(\\%%) & $\sigma$\\\\\n\hline\n" % (com, comm_size[com][0], comm_innerw[com] ) )
          for i in range(20):
            f_out.write(" &  & \\\\\n")
        #TK
        f_out.write("\hline\n\hline\nTitle Words & f(\\%) & $\sigma$\\\\\n\hline\n")
        if com in stuffTK:
          for i in range(len(stuffTK[com])):
            f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffTK[com][i][0], stuffTK[com][i][1], stuffTK[com][i][2]) )
          for i in range(max(0,10-len(stuffTK[com]))):
            f_out.write(" &  & \\\\\n")
        else:
          for i in range(10): f_out.write(" &  & \\\\\n")
        #J
        f_out.write("\hline\n\hline\nJournal & f(\\%) & $\sigma$\\\\\n\hline\n")
        if com in stuffJ:
          for i in range(len(stuffJ[com])):
            f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffJ[com][i][0][0:min(30,len(stuffJ[com][i][0])-1)], stuffJ[com][i][1], stuffJ[com][i][2]) )
          for i in range(max(0,10-len(stuffJ[com]))):
            f_out.write(" &  & \\\\\n")
        else:
          for i in range(10): f_out.write(" &  & \\\\\n")
        f_out.write("\hline\n\end{tabular}\n}\n")
        #f_out.write("\hline\n\end{tabular}\n}\n\end{table}\n\n")
        #I
        f_out.write("{\scriptsize\\begin{tabular}{|l r r|}\n\hline\nInstitution & f(\\%) & $\sigma$\\\\\n\hline\n")
        if com in stuffI:
          for i in range(len(stuffI[com])):
            if len(stuffI[com][i][0]) < 30:
              f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffI[com][i][0], stuffI[com][i][1], stuffI[com][i][2]) )
            else:
              aux = stuffI[com][i][0].rfind(' ')
              while aux > 30: 
                 aux = stuffI[com][i][0][0:aux].rfind(' ')
              f_out.write("%s &  & \\\\\n" % ( stuffI[com][i][0][0:aux] ) )
              f_out.write("$\quad$%s & %1.2f & %1.2f\\\\\n" % ( stuffI[com][i][0][aux:], stuffI[com][i][1], stuffI[com][i][2]) )
          for i in range(max(0,20-len(stuffI[com]))):
            f_out.write(" &  & \\\\\n")
        else:
          for i in range(20): f_out.write(" &  & \\\\\n")
        #C
        f_out.write("\hline\n\hline\nCountry & f(\\%) & $\sigma$\\\\\n\hline\n")
        if com in stuffC:
          for i in range(len(stuffC[com])):
            f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffC[com][i][0], stuffC[com][i][1], stuffC[com][i][2]) )
          for i in range(max(0,10-len(stuffC[com]))):
            f_out.write(" &  & \\\\\n")
        else:
          for i in range(10): f_out.write(" &  & \\\\\n")
        #A
        f_out.write("\hline\n\hline\nAuthor & f(\\%) & $\sigma$\\\\\n\hline\n")
    
        if com in stuffA:
          #@MM dump first authors (already sorted)
          #print "AUT; ", com, ";", stuffA[com][0][0], ";", stuffA[com][0][1], ";", stuffA[com][0][2]
          #aut[com]=(stuffA[com][0][0], stuffA[com][0][1], stuffA[com][0][2])
    
          for i in range(len(stuffA[com])):
            f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffA[com][i][0], stuffA[com][i][1], stuffA[com][i][2]) )
          for i in range(max(0,10-len(stuffA[com]))):
            f_out.write(" &  & \\\\\n")
        else:
          for i in range(10): f_out.write(" &  & \\\\\n")
        f_out.write("\hline\n\end{tabular}\n}\n")
        #R
        f_out.write("{\scriptsize\\begin{tabular}{|l r r|}\n\hline\nReference & f(\\%) & $\sigma$\\\\\n\hline\n")
        if com in stuffR:
          for i in range(len(stuffR[com])):
            if len(stuffR[com][i][0]) < 50:
              f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffR[com][i][0], stuffR[com][i][1], stuffR[com][i][2]) )
            elif len(stuffR[com][i][0]) < 90:
              aux = stuffR[com][i][0].rfind(' ')
              while aux > 50: 
                 aux = stuffR[com][i][0][0:aux].rfind(' ')
              f_out.write("%s &  & \\\\\n" % ( stuffR[com][i][0][0:aux] ) )
              f_out.write("$\quad$%s & %1.2f & %1.2f\\\\\n" % ( stuffR[com][i][0][aux:], stuffR[com][i][1], stuffR[com][i][2]) )
            else:
              aux1 = stuffR[com][i][0].rfind(' ')
              while aux1 > 90: 
                 aux1 = stuffR[com][i][0][0:aux1].rfind(' ')
              aux2 = stuffR[com][i][0][0:aux1].rfind(' ')
              while aux2 > 50: 
                 aux2 = stuffR[com][i][0][0:aux2].rfind(' ')
              f_out.write("%s &  & \\\\\n" % ( stuffR[com][i][0][0:aux2] ) )
              f_out.write("$\quad$%s &  & \\\\\n" % ( stuffR[com][i][0][aux2:aux1] ) )
              f_out.write("$\quad$%s & %1.2f & %1.2f\\\\\n" % ( stuffR[com][i][0][aux1:], stuffR[com][i][1], stuffR[com][i][2]) )
          for i in range(max(0,20-len(stuffR[com]))):
            f_out.write(" &  & \\\\\n")
        else:
          for i in range(20): f_out.write(" &  & \\\\\n")
        #RJ
        f_out.write("\hline\n\hline\nRefJournal & f(\\%) & $\sigma$\\\\\n\hline\n")
        if com in stuffRJ:
          for i in range(len(stuffRJ[com])):
            if len(stuffRJ[com][i][0]) < 50:
              f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffRJ[com][i][0], stuffRJ[com][i][1], stuffRJ[com][i][2]) )
            else:
              aux = stuffRJ[com][i][0].rfind(' ')
              while aux > 50: 
                 aux = stuffRJ[com][i][0][0:aux].rfind(' ')
              f_out.write("%s &  & \\\\\n" % ( stuffRJ[com][i][0][0:aux] ) )
              f_out.write("$\quad$%s &  %1.2f & %1.2f \\\\\n" % ( stuffRJ[com][i][0][aux:], stuffRJ[com][i][1], stuffRJ[com][i][2] ) )
          for i in range(max(0,10-len(stuffRJ[com]))):
            f_out.write(" &  & \\\\\n")
        else:
          for i in range(10): f_out.write(" &  & \\\\\n")
        #S
        if com in stuffS:
            f_out.write("\hline\n\hline\nSubject & f(\\%) & $\sigma$\\\\\n\hline\n")
            for i in range(len(stuffS[com])):
              f_out.write("%s & %1.2f & %1.2f\\\\\n" % ( stuffS[com][i][0], stuffS[com][i][1], stuffS[com][i][2]) )
            for i in range(max(0,10-len(stuffS[com]))):
              f_out.write(" &  & \\\\\n")
            #else:
            #  for i in range(10): f_out.write(" &  & \\\\\n")
        ###
        f_out.write("\hline\n\end{tabular}\n}\n\end{table}\n\n")
        
    #.. end
    f_out.write("\end{landscape}\n\n\end{document}\n")
    f_out.close()
    if verbose: print"..Communities characteristics extracted in .tex 'IDCards' file"
    
    #return
    
def cmpval(x,y):
    if (x[1][1]<y[1][1]) or (x[1][1]==y[1][1] and x[1][0]>y[1][0]):
        return -1
    elif x[1][1]==y[1][1] and x[1][0]==y[1][0]:
        return 0    
    else:
        return 1

def cmpvalb(x,y):
    if (x[1][2]<y[1][2]) or (x[1][2]==y[1][2] and x[1][0]>y[1][0]):
        return -1
    elif x[1][2]==y[1][2] and x[1][0]==y[1][0]:
        return 0    
    else:
        return 1      
        
def extract_labels(stuff):
    comm_label=dict()
    # label 0 is the most freq
    comm_label[0] = stuff[0][0]
    # label 1 is the most sign
    f=0
    for ff in range(len(stuff)):
      if stuff[ff][2]>stuff[f][2]: 
        f=ff 
    comm_label[1] = stuff[f][0]
    # label 2 has the max freq * sign
    f=0
    for ff in range(len(stuff)):
      if stuff[ff][1]*stuff[ff][2]>stuff[f][1]*stuff[f][2]: 
        f=ff 
    comm_label[2] = stuff[f][0]
    # label 3 is the most freq with sign > 1
    fff=0; 
    while (fff < len(stuff)-1) and (stuff[fff][2] <=0):
      fff+=1 
    while (fff < len(stuff)-1) and (stuff[fff][2] <=1):
      fff+=1   
    comm_label[3] = stuff[fff][0]
    #
    return comm_label

if __name__ == "__main__":
    main()
