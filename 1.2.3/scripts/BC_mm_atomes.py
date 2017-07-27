#! /usr/bin/env python

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2012
   All rights reserved.
   BSD license.
   .... If you are using these scripts, please cite our "Scientometrics" paper:
   .... S Grauwin, P Jensen, Mapping Scientific Institutions. Scientometrics 89(3), 943-954 (2011)
"""

# usage: BC.py -i DIR [-o DIR] [-p] [-v]
# if the output dir is not specified, the code use the input folder as an output folder 
# use option -p when the partitions have already been computed and you want to use them
# use option -v for verbose informations

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
import community_atoms
from subprocess import call, Popen, PIPE


#AFAF, usare media*sigma per autori e keyword
#      output gexf
## ##################################################
## ##################################################
## ##################################################
def BC_network(in_dir,out_dir,verbose,presaved,ids_fn,aut_fn):
  ## INI
  t1=time.time()

  ## INPUT DATA
  if verbose: print "..Initialize"

  src1  = os.path.join(in_dir, "articles.dat") 
  src6  = os.path.join(in_dir, "countries.dat") 
  src5  = os.path.join(in_dir, "references.dat")
  
  ## lire les infos articles
  Ymin = 2100; Ymax = 1900 # store the min and max publication year 
  nR = dict() # store the number of refs of the articles
  firstCountry = dict()
  pl = Utils.Article()
  pl.read_file(src1)  
  nb_art = len(pl.articles) # store the number of articles within database
  for l in pl.articles:
      firstCountry[l.id]=""
      nR[l.id] = 0
      if (l.year > 1900 and l.year < 2100): 
          if(l.year > Ymax): Ymax=l.year
          if(l.year < Ymin): Ymin=l.year

  ##lire premier pays de chaque article
  ql = Utils.Country()
  ql.read_file(src6)
  for l in ql.countries:
    if l.rank==0: firstCountry[l.id]=l.country  

  ## CREATE BC WEIGHT TABLE
  if verbose: print "..Create the 'Bibliographic Coupling' weight table"
 
  ref_table = dict() # store the id of articles using a given ref
  BC_table = dict() # store the number of common refs between pairs of articles

  if verbose: print "....loading refs table"
  pl = Utils.Ref()
  pl.read_file(src5)  
  for l in pl.refs:
      foo = l.firstAU + ', ' + str(l.year) + ', ' + l.journal + ', ' + l.volume + ', ' + l.page 
      if foo in ref_table: ref_table[foo].append( l.id )
      else: ref_table[foo] = [l.id]
      nR[l.id] += 1

  if verbose: print "....detecting common references"
  for foo in ref_table:
    if len(ref_table[foo]) > 1:
      for i in ref_table[foo]:
          for j in ref_table[foo]:
              if i<j:
                  if i not in BC_table: BC_table[i] = dict()
                  if j not in BC_table[i]: BC_table[i][j] = 0
                  BC_table[i][j] += 1 

  t2=time.time()
  print '..time needed until now: %ds' % (t2-t1)

  # choose threshold
  #confirm = 'n'; thr=2;
  confirm = 'n'; thr=2;
  
  while confirm != 'y':
    if thr == 1: print "Keep BC links between articles sharing at least %d reference" % (thr)
    else: print "Keep BC links between articles sharing at least %d references" % (thr)
    #confirm = raw_input("Confirm (y/n): ")
    confirm = 'y'
    while confirm not in ['n','y']:
      confirm = raw_input("...typing error!\n Confirm (y/n): ")
    if confirm == 'n':
      thr = input("threshold for BC links -- articles should share at least ? references:")

  bcthr = thr

  ##############################
  ## BC COMMUNITIES
  if verbose: print "..BC communities"
  #... define BC network
  if verbose: print "....define graph in networkx format"
  G=nx.Graph()
  TOTW=0
  for i in BC_table:
    for j in BC_table[i]:
      if BC_table[i][j] >= bcthr:
        w_ij = (1.0 * BC_table[i][j]) / math.sqrt(nR[i] * nR[j])
        TOTW+=w_ij
        G.add_edge(i, j, weight=w_ij)

  # compute stuff
  h=nx.degree(G).values()
  avg_degree=sum(h)*1.0/len(h)
  avg_weight=2*TOTW*1.0/(len(G.nodes())*(len(G.nodes())-1))

  #hw = nx.degree(G, weight=).values()
  #avg_w_degree = sum(h)
  density = nx.density(G) 

  #... output infos
  print "....There are %d articles in the database" % (nb_art)
  print "....There are %d articles in the BC network\n......(ie sharing at least %d reference(s) with another article)" % (len(G.nodes()), bcthr)
  #print "....BC networks: Average degree: %.3f, average weight: 1/%d" % (avg_degree, round(1/avg_weight))
  #MM
  print "....BC networks: Average degree: %.3f, average weight: 1/%d, density: %.3f" % (avg_degree, round(1/avg_weight), density)

  t2=time.time()
  print '..time needed until now: %ds' % (t2-t1)

  ##******************************************************************
  #algo_method = 'python' c++ thrown away
  ##******************************************************************

  if True: #was used to test for algo_method
    if verbose: print "....computing communities with networkx Louvain algo"
    dendogram = community_atoms.generate_dendogram(G, part_init=None)

    #... second louvain partition
    """
    if verbose: print "....computing level %d sub-communities with Louvain algo" % (len(dendogram))
    part = community_atoms.partition_at_level(dendogram,len(dendogram)-1)
    part2=part.copy()
    for com in set(part.values()) :
      list_nodes = [nodes for nodes in part.keys() if part[nodes] == com]
      # split communities of size > 100 
      if len(list_nodes) > 100: 
        H = G.subgraph(list_nodes)
        dendo2 = community_atoms.generate_dendogram(H, part_init=None)
        partfoo = community_atoms.partition_at_level(dendo2,len(dendo2)-1)

         # add prefix code
        for aaa in partfoo.keys():
          partfoo[aaa] = (com+1)*1000 + partfoo[aaa] ###UGH.
        mod = community_atoms.modularity(partfoo, H)
        nb_comm = len(set(partfoo.values()))
        if verbose: print "..... splitting comm %d (N=%d records) in %d sub-comm, Q=%.3f" % (com, len(list_nodes), nb_comm, mod)

        #MM hierarchical test: if (i.): mod > 0.4 ; then if (ii.) subcommunities avg size > 10 ##
        subcomm_q = mod
        
        print "....... hierarchical test: subcomm Q = %.3f" % (subcomm_q)
        if (subcomm_q <0.4):
          print "......... First test (q>.4) failed. Not splitting."
        else: 
          print "......... First test (q>.4) passed."

          print "average_comm_size: "
          nb_subcomm = {}

          for v in partfoo.itervalues():
            if not (nb_subcomm.has_key(v)):
              nb_subcomm[v]=0
            else:
              nb_subcomm[v]+=1

          print "*****nb_subcomm", nb_subcomm
          average_comm_size=len(list_nodes)*1.0/nb_comm
          print "AVERAGE COMM SIZE: %.3f", average_comm_size

          if (average_comm_size > 10.0):
            print "......... Second test (avgsize>10) passed. Splitting."
            part2.update(partfoo)
          else:
            print "......... Second test (avgsize>10) failed. Not splitting."
    """


    #... save partitions
    fooname="partitions_thr%d.txt" % (bcthr)
    filename = os.path.join(out_dir, fooname)
    f_out = open(filename,"w")
    louvain_partition=dict()
    for lev in range(len(dendogram)):
      louvain_partition[lev]=community_atoms.partition_at_level(dendogram, lev) 
    #.. I want communtity labels starting from 1 instead of 0 for top level  
    for k in louvain_partition[len(dendogram)-1].keys():
      louvain_partition[len(dendogram)-1][k] += 1
    #louvain_partition[len(dendogram)]=part2
    f_out.write('%s' % json.dumps(louvain_partition))
    f_out.close()

    t2=time.time()
    print '..time needed until now: %ds' % (t2-t1)
    

  else:
    if verbose: print "....uploading previously computed communities"
    #... upload partition
    fooname="partitions_thr%d.txt" % (bcthr)
    filename = os.path.join(out_dir, fooname)
    if not os.path.isfile(filename):
      print '....file %s does not exists' % filename
      return
    ffin = open(filename,'r')
    foo = ffin.read()
    lines = foo.split('\n')
    #for line in lines:
    louvain_partition = json.loads(lines[0])
    num_levels=len(louvain_partition)-1
    print "....%d+1 levels" % num_levels
    #... convert keys back into integer (json dump put them in strings)
    auxlist = louvain_partition.keys()
    for k in auxlist:
      louvain_partition[int(k)]=louvain_partition[k] 
      del louvain_partition[k]
      auxlistb = louvain_partition[int(k)].keys()   
      for kk in auxlistb:
        louvain_partition[int(k)][int(kk)]=louvain_partition[int(k)][kk] 
        del louvain_partition[int(k)][kk]  

  #... output infos
  num_levels=len(louvain_partition)-1 ####MM TOP, no SUBTOP
  for level in range(len(louvain_partition)):
    if level < len(louvain_partition)-1:
      txt = "level %d" % (level + 1)
    else:
      txt = "level %d SUB" % (level)
    part = louvain_partition[level]  
    mod = community_atoms.modularity(part, G)
    nb_comm = len(set(part.values()))
    size_sup10 = 0; size_sup100 = 0; size_sup5 = 0; #communities_caracteristics(partition, thr, level)
    for com in set(part.values()) :
      list_nodes = [nodes for nodes in part.keys() if part[nodes] == com]
      if len(list_nodes) > 100: size_sup100 += 1
      if len(list_nodes) > 10: size_sup10 += 1
      #MM
      if len(list_nodes) > 5: size_sup5 += 1
    print "....%s: %d communities [%d with size > 5, %d with size > 10, %d with size > 100], modularity Q=%1.6f" % (txt, nb_comm, size_sup5, size_sup10, size_sup100, mod)

  #MM
  #hlvl = len(louvain_partition)-2 #non-sub
  #hlvl_louvpart = louvain_partition[hlvl]
  #hlvlq = community.modularity(hlvl_louvpart, G)
  #hlvlncomm = len(set(hlvl_louvpart.values()))
  #print "Highest non-sub level, Modularity, #comm, comm sized>5: %d ; %1.6f ; %d ; %d" % (hlvl+1, hlvlq, hlvlncomm, size_sup5)

  hlvl_s = len(louvain_partition)-1 #sub
  hlvl_louvpart_s = louvain_partition[hlvl_s]
  hlvlq_s = community_atoms.modularity(hlvl_louvpart_s, G)
  hlvlncomm_s = len(set(hlvl_louvpart_s.values()))
  print "Highest     sub level, Modularity, #comm, comm sized>5: %d ; %1.6f ; %d ; %d" % (hlvl_s+1, hlvlq_s, hlvlncomm_s, size_sup5)

  t2=time.time()
  print '..time needed until now: %ds' % (t2-t1)

  ############################################################
  ############################################################
  ## EXTRACT COMMUNITIES?
  #confirm = raw_input("..Do you want to create a gephi file and an 'IDcards' document at a communities level? (y/n): ")
  confirm = 'y'
  if confirm == 'y':
    ##############################
    ## WHICH EXTRACTION ?
    print "..BC communities extraction"
    #
    #level = input("......level you want to extract:")
    #thr  = input("......keep communities of size > to:")
    #LEVEL LOUVAIN SET HERE @#@#@#
    thr = 5
    level = num_levels;
    confirm = 'n'; 
    #level = num_levels; thr = 50;
    
    while confirm != 'y':
      part=louvain_partition[(level - 1)].copy()
      txt = 'level%d' % level
      textxt = 'level %d/%d' % (level,num_levels)
      if level == num_levels:
        textxt = 'top level'
        #yesno = raw_input("......split large (>100) clusters? (y/n): ")
        yesno = 'n' ###want to compute subtop partitions?
        if yesno == 'y':
          part2=louvain_partition[num_levels] # this is the all split communities 
          txt = 'level%d_sub' % level
          textxt = 'sub-top level'
          Qthr = input("......keep split with Q >:")
          Nthr = input("......keep split with N >:")
          for com in set(part.values()):
            list_nodes = [nodes for nodes in part.keys() if part[nodes] == com]
            # split communities of size > 100 
            if len(list_nodes) > 100: 
              H = G.subgraph(list_nodes)
              partfoo = dict((k, part2[k]) for k in list_nodes)
              mod = community_atoms.modularity(partfoo, H)
              nb_comm = len(set(partfoo.values()))
              if (mod >= Qthr) and (len(list_nodes) >= Nthr):
                if verbose: print "..... community %d is split in %d sub-communities, Q=%.3f" % (com, nb_comm, mod)
                part.update(partfoo)
      #..infos          
      nb_comm = len(set(part.values()))
      size_sup_thr = 0; n_sup_thr = 0;
      for com in set(part.values()) :
        list_nodes = [nodes for nodes in part.keys() if part[nodes] == com]
        if len(list_nodes) > thr: 
          n_sup_thr += len(list_nodes)
          size_sup_thr += 1
      print "....Extraction of %s BC communities with size > %d\n......(%d articles gathered in %d communities):" % (txt, thr, n_sup_thr, size_sup_thr)
      #confirm = raw_input("....do you confirm? (y/n): ")
      confirm = 'y'
      # ..confirm
      if confirm == 'n':
        level = input("......level you want to extract:")
        thr  = input("......keep communities of size > to:")

    which ='x'
    while (which != 's' and which !='d'):  
      #which = raw_input("......in the document presenting the clusters characteristics, rank the clusters by descending size or density? (s/d):")
      which = 's'

    ##############################
    ## DIRECTION + GLUING REFERENCES
    if verbose: print"..computing direction and gluing references"
    # ini
    fooname="glue_thr%d_%s.txt" % (bcthr,txt)
    filename = os.path.join(out_dir, fooname)
    f_out = open(filename,"w") 
    # do
    [stuff, direction] = BCUtils.direction_and_glue(G,in_dir,part,thr,verbose,txt)
    # output
    if verbose: print"..Output most gluing references"
    for i in range(len(stuff)):
      f_out.write('%s\t%.3f\n' % (stuff[i][0],stuff[i][1]))
    # end
    f_out.close() 
    if verbose: print"..Done!\n"
    t2=time.time()
    print '..time needed until now: %ds' % (t2-t1)  

    ##############################
    ## COMMUNITIES CARACTERISTICS
    if verbose: print"..Computing communities caracteristics"

    #.. ini
    fooname="BCcomm_ID_Cards_thr%d_%s.tex" % (bcthr,txt)
    filename = os.path.join(out_dir, fooname)
    f_out = open(filename,"w")
    f_out.write("\documentclass[a4paper,11pt]{report}\n\usepackage[english]{babel}\n\usepackage[latin1]{inputenc}\n\usepackage{amsfonts,amssymb,amsmath}\n\usepackage{pdflscape}\n\usepackage{color}\n\usepackage{graphicx}\n\n\\addtolength{\evensidemargin}{-60pt}\n\\addtolength{\oddsidemargin}{-60pt}\n\\addtolength{\\textheight}{80pt}\n\n\\title{{\\bf Communities ID Cards}}\n\date{\\begin{flushleft}This document gathers the ``ID Cards'' of the BC communities found within the studied database.\\\\\n The BC network was built by linking pairs of articles sharing at least %d references - %d out %d articles are in the network. The %d communities presented here correspond to the ones found in the %s level grouping at least %d articles. They correspond to a total of %d articles. \\\\\n These ID cards displays the most frequent keywords, subject categories, journals of publication, institution, countries, authors, references and reference journals of the articles of each community. The significance of an item $\sigma = \sqrt{N} (f - p) / \sqrt{p(1-p)}$ [where $N$ is the number of articles within the community and $f$ and $p$ are the proportion of articles respectively within the community and within the database displaying that item ] is also given (for example $\sigma > 5$ is really highly significant). Each community is labelled by the most frequent keyword with positive $\sigma$.\\\\\n\\vspace{1cm}\n\copyright Sebastian Grauwin - BiblioTools3.0 (September 2014) \end{flushleft}}\n\n\\begin{document}\n\\begin{landscape}\n\maketitle\n" % (bcthr, len(G.nodes()), nb_art, size_sup_thr, textxt, thr, n_sup_thr))

    f_out.write("\n\n\clearpage\\begin{figure}[h!]\n\\begin{center}\n%%\includegraphics[width=1.3\\textwidth]{xxx}\n\caption{{\\bf Top communities network.} The sizes of the nodes correspond to the number of articles in each community. Their colors correspond to the top community they belong to and the labels correspond to those reported in the tables below. The edges reflects shared references between communities. The thicker an edge between two communities, the more references they share. The arrows reflect a `who feed who' relationship: an arrow pointing from A to B indicates that the references shared by A and B are more similar to the references from A.}\n\end{center}\n\end{figure}\n\n")

    #... partition
    partition = part
    list_nodes = dict();
    for com in set(partition.values()) :
      list_nodes[com] = [nodes for nodes in partition.keys() if partition[nodes] == com]

	#json dictionaries for labels (main kw, authors)
    ids={}
    aut={}


    #... about top partition
    top_part = louvain_partition[len(louvain_partition)-2]
    # com to top comm info
    top_com = dict()
    for com in list_nodes:
      top_com[com]=top_part[list_nodes[com][0]]
    # inner weight of top com
    list_top = dict();
    topcomm_innerw = dict(); topcomm_size = dict();
    for topcom in set(top_part.values()) :
      list_top[topcom] = [nodes for nodes in top_part.keys() if top_part[nodes] == topcom]
      size = len(list_top[topcom])
      W = 0;
      for id1 in list_top[topcom]:
        for id2 in list_top[topcom]:
          if id2 > id1 and id2 in G.edge[id1]: 
            W += G.edge[id1][id2]['weight']
      W *= 2.0 / (size * (size -1))
      topcomm_innerw[topcom]= 1.0 / W
      topcomm_size[topcom] = [size,1]


    #... quantitative stuff about partition
    comm_innerw = dict(); comm_size = dict(); 
    for com in list_nodes:
      size = len(list_nodes[com])
      W = 0;
      for id1 in list_nodes[com]:
        for id2 in list_nodes[com]:
          if id2 > id1 and id2 in G.edge[id1]: 
            W += G.edge[id1][id2]['weight']
      W *= 2.0 / (size * (size -1))
      comm_innerw[com]= 1.0 / W
      comm_size[com] = [size,topcomm_innerw[top_com[com]],topcomm_size[top_com[com]]]
    #...order communities by innerw (or size, depending on previously defined parameter "which") of their top and then by size 
    Lcomm_size = comm_size.items()
    if which =='s': Lcomm_size.sort(cmpval)
    if which =='d': Lcomm_size.sort(cmpvalb)

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
          ids[com]=(comm_label[com][0], comm_label[com][3], author_mostfreq[0], comm_size[com][0])

          f_out.write("\clearpage\n\n\\begin{table}[!ht]\n\caption{The community ``%s'' (id %d) contains $N = %d$ articles. Its average internal link weight is $<\omega_{in}> \simeq 1/%d$ }\n\\textcolor{white}{aa}\\\\\n{\scriptsize\\begin{tabular}{|l r r|}\n\hline\nKeyword & f(\\%%) & $\sigma$\\\\\n\hline\n" % (comm_label[com][3], com, comm_size[com][0], comm_innerw[com] ) )
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
          aut[com]=(stuffA[com][0][0], stuffA[com][0][1], stuffA[com][0][2])

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
    if verbose: print"..Communities caracteristics extracted in .tex 'IDCards' file"

	## DUMP JSON FILES
    #json.dump(ids, open(ids_fn,'wb'))
    #json.dump(aut, open(aut_fn,'wb'))

    ##############################
    ## OUTPUT GEPHI FILES

    #... output gephi
    if verbose: print "..Preparing gephi gdf file for BC communities network"

    ## ... ini
    name = "BC_thr%d_comm_%s.gdf" % (bcthr, txt)
    dst = os.path.join(out_dir, name)

    Gpart = nx.Graph()

    f_gephi = open(dst,'w')
    ## ... prep nodes
    if verbose: print "....nodes"
    f_gephi.write("nodedef>name VARCHAR,label_f VARCHAR,label_s VARCHAR,label_fs VARCHAR,label VARCHAR,topcom VARCHAR,colortop VARCHAR,size DOUBLE,inv_innerweight DOUBLE\n")
    for com in comm_size:
      if comm_size[com][0] > thr: 
          
          f_gephi.write("%d,'%s','%s','%s','%s',%d,color%d,%d,%1.0f\n" % (com, comm_label[com][0],comm_label[com][1],comm_label[com][2],comm_label[com][3],top_com[com],top_com[com],comm_size[com][0], comm_innerw[com]) )    

          Gpart.add_node(int(com), {'comm_label0': comm_label[com][0],\
                             'comm_label1': comm_label[com][1],\
                             'comm_label2': comm_label[com][2],\
                             'comm_label3': comm_label[com][3],\
                             'main_author': aut[com][0],\
                             'top_com': top_com[com],\
                             'comm_size': comm_size[com][0],\
                             'comm_innerw': comm_innerw[com]})

    ## ... prep edges
    if verbose: print "....edges"
    f_gephi.write("edgedef>node1 VARCHAR,node2 VARCHAR,num_links DOUBLE,weight DOUBLE,logweight DOUBLE,cos1 DOUBLE,cos2 DOUBLE\n")
    for com1 in list_nodes:
      for com2 in list_nodes:
        size1 = len(list_nodes[com1]); size2 = len(list_nodes[com2]);
        if size1 > thr and size2 > thr and com1 > com2:
          W = 0; N = 0;
          for id1 in list_nodes[com1]:
            for id2 in list_nodes[com2]:
              if id2 in G.edge[id1]: 
                N += 1
                W += G.edge[id1][id2]['weight']
          W *= 1000.0 / (size1 * size2)
          cos1=direction[(com1,com2)][0]
          cos2=direction[(com1,com2)][1]          
          if W > 0.000001:
            if cos1 < cos2:
              f_gephi.write("%d,%d,%d,%1.9f,%1.2f,%.4f,%.4f\n" % (com1, com2, N, W, 6 + math.log(W)/math.log(10),cos1,cos2))
            else:
              f_gephi.write("%d,%d,%d,%1.9f,%1.2f,%.4f,%.4f\n" % (com2, com1, N, W, 6 + math.log(W)/math.log(10),cos2,cos1))  
              
            Gpart.add_edge(int(com1), int(com2), {'N': N, 'weight': W,\
                                        'logweight': 6 + math.log(W)/math.log(10),\
                                        'cos1': cos1, 'cos2': cos2})
    ## ... end
    f_gephi.close() 

    #DUMP gexf file (all nodes)
    #make an enriched graph
    gexf_part_fileout = os.path.join(out_dir, 'partitions.gexf')
    nx.write_gexf(Gpart, gexf_part_fileout)


    if verbose: print"..Done!\n"
    # ##
    # ##

  ############################################################
  ############################################################
  ## OUTPUT A GEPHI FILE AT THE NODE LEVEL
  ##... output the BC networks?  
  #confirm = raw_input("..Do you want to create a gephi file at the article level? (y/n): ") 
  confirm = 'y'
  if confirm == 'y':
    while confirm == 'y':
      p=louvain_partition[num_levels-1]
      #confirm2 = raw_input("....do you want to keep the whole network, ie %d nodes / %d edges? (y/n): " % (len(G.nodes()), len(G.edges()) ))
      confirm2 = 'y'
      if confirm2 =='y':
        nodes_to_keep=[n for n in p]
        name = "network_all.gdf"
      else:
        commtokeep=input("....keep only top community number:")
        nodes_to_keep=[n for n in p if p[n]==commtokeep]
        name = "network_topcomm%d.gdf" % commtokeep
      ## ... ini
      dst = os.path.join(out_dir, name)
      f_gephi = open(dst,'w')
      ## ... prep nodes
      if verbose: print "....nodes"
      f_gephi.write("nodedef>name VARCHAR,label VARCHAR,top_com VARCHAR,subtop_com VARCHAR,bottom_com VARCHAR,firstAU VARCHAR,journal VARCHAR,year VARCHAR,nb_refs DOUBLE,citations DOUBLE,country VARCHAR\n")
      ##appelle liste articles, class Article format specifique dans Utils.py
      ##src1 = articles.dat, chaque article ==l pointeur   

      pl = Utils.Article()
      pl.read_file(src1)
      for l in pl.articles:
        if l.id in nodes_to_keep: 
          topcom = str(louvain_partition[num_levels-1][l.id])
          subtopcom = str(louvain_partition[num_levels][l.id])
          bottomcom = str(louvain_partition[0][l.id])
          foo = l.firstAU + ', ' + l.journal + ', ' + str(l.year)
          f_gephi.write("%d,'%s',%s,%s,%s,'%s','%s',%d,%d,%s,'%s'\n" % (l.id, foo, topcom, subtopcom, bottomcom, l.firstAU, l.journal,l.year,nR[l.id],l.times_cited,firstCountry[l.id]) )  

          (u_id, u_descr, u_part, u_auth, u_journ, u_year, u_nr, u_cit, u_country) = \
              (l.id, foo, topcom, l.firstAU, l.journal,l.year,nR[l.id],l.times_cited,firstCountry[l.id]) 
              
          G.node[u_id]['part']=u_part
          G.node[u_id]['descr']=u_descr
          G.node[u_id]['part']=u_part
          G.node[u_id]['auth']=u_auth
          G.node[u_id]['journ']=u_journ
          G.node[u_id]['year']=u_year
          G.node[u_id]['nr']=u_nr
          G.node[u_id]['cit']=u_cit
          G.node[u_id]['country']=u_country


      ## ... prep edges
      if verbose: print "....edges"
      f_gephi.write("edgedef>node1 VARCHAR,node2 VARCHAR,weight DOUBLE,nb_comm_refs DOUBLE")
      H = G.subgraph(nodes_to_keep)
      for e in H.edges(data=True):
        w_ij = e[2]['weight']
        i=min(e[0],e[1])
        j=max(e[0],e[1])
        f_gephi.write("\n%d,%d,%f,%d" % (i, j, w_ij, BC_table[i][j]))
        
        G.add_edge(i, j, {'weight': w_ij, 'nb_comm_refs': BC_table[i][j]})
      ## ... end
      f_gephi.close()
      
      #DUMP gexf file (all nodes)
      #make an enriched graph
      gexf_fileout = os.path.join(out_dir, 'nodes.gexf')
      nx.write_gexf(G, gexf_fileout)
      
      #confirm = raw_input("..Do you want to create another gephi file at the article level? (y/n): ") 
      confirm = 'n'

  ##
  ##    
  if verbose: print"..Done!\n"


  ## ###################################
  ## END
  t2=time.time()
  print 'total time needed: %ds' % (t2-t1)
  return

## ##################################################
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

## ##################################################

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
## ##################################################
## ##################################################

def main():
# usage: BC.py [-h] [--version] -i DIR -o DIR[-v]
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --version             show program's version number and exit
#   -i DIR, --input_dir DIR input directory name
#   -g 
#   -o DIR, --output_dir DIR input directory name
  # Parse line options.
  # Try to have always the same input options
  parser = argparse.ArgumentParser(description = 'parser')

  parser.add_argument('--version', action='version', version='%(prog)s 1.1')
  
  parser.add_argument("-i", "--input_dir", nargs=1, required=True,
          action = "store", dest="in_dir",
          help="input directory name",
          metavar='DIR')

  parser.add_argument("-o", "--output_dir", nargs=1, required=False,
          action = "store", dest="out_dir", 
          default = 'Desktop/',
          help="output directory name",
          metavar='DIR')

  parser.add_argument("-p", "--presaved",
          action = "store_true", dest="presaved",
          default = False,
          help="presaved mode [default %(default)s]")
          
  parser.add_argument("-v", "--verbose",
          action = "store_true", dest="verbose",
          default = False,
          help="verbose mode [default %(default)s]")

  #Analysis of input parameters
  args = parser.parse_args()
  
  if not os.path.exists(args.in_dir[0]):
      print "Error: Input directory does not exist: ", args.in_dir[0]
      exit()

  if args.out_dir == 'Desktop/':
    args.out_dir = args.in_dir

  if not os.path.exists(args.out_dir[0]):
      print "Error: Output directory does not exist: ", args.out_dir[0] 
      exit()

  print args.in_dir
  ids_fn = 'ids'+str(args.in_dir[0])+'.json'
  aut_fn = 'aut'+str(args.in_dir[0])+'.json'

  ##      
  BC_network(args.in_dir[0],args.out_dir[0],args.verbose,args.presaved, ids_fn, aut_fn)

  return


    
## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################
