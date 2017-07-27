#! /usr/bin/env python

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2012
   All rights reserved.
   BSD license.
   .... If you are using these scripts, please cite our "Scientometrics" paper:
   .... S Grauwin, P Jensen, Mapping Scientific Institutions. Scientometrics 89(3), 943-954 (2011)
"""

import os
import sys
import glob
import numpy
import argparse
import math
import Utils

## ##################################################
## ##################################################
## ##################################################
def comm_tables(in_dir,partition,thr,verbose):

  ## INPUT DATA
  src1  = os.path.join(in_dir, "articles.dat") 
  src2  = os.path.join(in_dir, "authors.dat")
  src3  = os.path.join(in_dir, "keywords.dat")
  src4  = os.path.join(in_dir, "subjects.dat")
  src5  = os.path.join(in_dir, "references.dat")
  src6  = os.path.join(in_dir, "countries.dat")
  src7  = os.path.join(in_dir, "institutions.dat")


  ## Communities sizes - we are mostly interested by articles within communities of size > thr
  comm_size = dict();
  for com in set(partition.values()) :
    list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]
    comm_size[com] = len(list_nodes)

  ## TREAT DATA
  #######
  # initialize
  list_comm = [] ## list of communities of size > thr
  list_id = dict(); ## list of id -- articles within the BC network and within a community of size > thr
  for com in comm_size: 
    if comm_size[com] > thr: list_comm.append( com )
  for id_art in partition:
    com = partition[id_art]
    if comm_size[com] > thr: list_id[id_art] = ''

  #######
  # KEYWORDS
  if verbose: print "....most frequent keywords"
  art_wK = dict(); ## lists article with keywords
  probaK = dict(); ## records the freq of each keyword in the whole database
  freqK  = dict(); ## "" in each community of size > thr
  for com in list_comm: freqK[com] = dict();

  # read data
  pl = Utils.Keyword()
  pl.read_file(src3)  
  for l in pl.keywords:
    #if (l.ktype == 'IK'): ###@MM - why IK only???
    if (l.ktype == 'IK' or l.ktype == 'AK'): ###@MM AK too (in CO2 1975 for e.g. no IK)
      # list the articles with keywords
      if l.id not in art_wK: art_wK[l.id] = ''
      # record the number of occurrence of a given K
      if l.keyword not in probaK: probaK[l.keyword] = 0
      probaK[l.keyword] += 1
      # record the occurrence of a given K within each community
      if l.id in list_id: 
        com = partition[l.id]
        if l.keyword not in freqK[com]: freqK[com][l.keyword] = 0
        freqK[com][l.keyword] += 1

  # extract the 20 most frequent keywords within each community, normalize their frequencies and compute their significance
  stuffK = dict();
  NK = len(art_wK);
  print "BCU - NK len=",NK
  if NK > 0:
    for com in freqK:
      cs = comm_size[com]
      stuffK[com] = dict()
      L = freqK[com].items()
      L.sort(cmpval)
      for i in range(min(20,len(L))):
        keyw = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaK[keyw] * 1.0 / NK 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffK[com][i] = [keyw, f*100, sigma]

  #######
  # TITLE WORDS
  if verbose: print "....most frequent title words"
  art_wTK = dict(); ## lists article with title words
  probaTK = dict(); ## records the freq of each word in the whole database
  freqTK  = dict(); ## "" in each community of size > thr
  for com in list_comm: freqTK[com] = dict();

  # read data
  pl = Utils.Keyword()
  pl.read_file(src3)  
  for l in pl.keywords:
    if (l.ktype == 'TK'):
      # list the articles with title words
      if l.id not in art_wTK: art_wTK[l.id] = ''
      # record the number of occurrence of a given TK
      if l.keyword not in probaTK: probaTK[l.keyword] = 0
      probaTK[l.keyword] += 1
      # record the occurrence of a given TK within each community
      if l.id in list_id: 
        com = partition[l.id]
        if l.keyword not in freqTK[com]: freqTK[com][l.keyword] = 0
        freqTK[com][l.keyword] += 1

  # extract the 10 most frequent title words within each community, normalize their frequencies and compute their significance
  stuffTK = dict();
  NTK = len(art_wTK);
  if NTK > 0:
    for com in freqTK:
      cs = comm_size[com]
      stuffTK[com] = dict()
      L = freqTK[com].items()
      L.sort(cmpval)
      for i in range(min(10,len(L))):
        tword = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaTK[tword] * 1.0 / NTK 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffTK[com][i] = [tword, f*100, sigma]

  #######
  # SUBJECTS
  if verbose: print "....most frequent subjects"
  art_wS = dict(); ## lists article with subjects
  probaS = dict(); ## records the freq of each subject in the whole database
  freqS  = dict(); ## "" in each community of size > thr
  for com in list_comm: freqS[com] = dict();

  # read data
  pl = Utils.Subject()
  pl.read_file(src4)  
  for l in pl.subjects:
      # list the articles with subjects
      if l.id not in art_wS: art_wS[l.id] = ''
      # record the number of occurrence of a given S
      if l.subject not in probaS: probaS[l.subject] = 0
      probaS[l.subject] += 1
      # record the occurrence of a given S within each community
      if l.id in list_id: 
        com = partition[l.id]
        if l.subject not in freqS[com]: freqS[com][l.subject] = 0
        freqS[com][l.subject] += 1
  #
  # extract the 10 most frequent subjects within each community, normalize their frequencies and compute their significance
  stuffS = dict();
  NS = len(art_wS);
  if NS > 0:
    for com in freqS:
      cs = comm_size[com]
      stuffS[com] = dict()
      L = freqS[com].items()
      L.sort(cmpval)
      for i in range(min(10,len(L))):
        subj = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaS[subj] * 1.0 / NS 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffS[com][i] = [subj.replace('&','\\&'), f*100, sigma]

  #######
  # JOURNALS
  if verbose: print "....most frequent journals"
  art_wJ = dict(); ## lists article with journals
  probaJ = dict(); ## records the freq of each journal in the whole database
  freqJ  = dict(); ## "" in each community of size > thr
  for com in list_comm: freqJ[com] = dict();

  # read data
  #
  pl = Utils.Article()
  pl.read_file(src1)  
  for l in pl.articles:
      # list the articles with journals
      if l.id not in art_wJ: art_wJ[l.id] = ''
      # record the number of occurrence of a given J
      if l.journal not in probaJ: probaJ[l.journal] = 0
      probaJ[l.journal] += 1
      # record the occurrence of a given J within each community
      if l.id in list_id: 
        com = partition[l.id]
        if l.journal not in freqJ[com]: freqJ[com][l.journal] = 0
        freqJ[com][l.journal] += 1
  #
  # extract the 10 most frequent journals within each community, normalize their frequencies and compute their significance
  stuffJ = dict();
  NJ = len(art_wJ);
  if NJ > 0:
    for com in freqJ:
      cs = comm_size[com]
      stuffJ[com] = dict()
      L = freqJ[com].items()
      L.sort(cmpval)
      for i in range(min(10,len(L))):
        jour = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaJ[jour] * 1.0 / NJ 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffJ[com][i] = [jour.replace('&','\\&'), f*100, sigma]

  #######
  # AUTHORS
  if verbose: print "....most frequent authors"
  art_wA = dict(); ## lists article with authors
  probaA = dict(); ## records the freq of each author in the whole database
  freqA  = dict(); ## "" in each community of size > thr
  for com in list_comm: freqA[com] = dict();

  # read data
  pl = Utils.Author()
  pl.read_file(src2)  
  for l in pl.authors:
      # list the articles with authors
      if l.id not in art_wA: art_wA[l.id] = ''
      # record the number of occurrence of a given A
      if l.author not in probaA: probaA[l.author] = 0
      probaA[l.author] += 1
      # record the occurrence of a given A within each community
      if l.id in list_id: 
        com = partition[l.id]
        if l.author not in freqA[com]: freqA[com][l.author] = 0
        freqA[com][l.author] += 1
  #
  # extract the 10 most frequent authors within each community, normalize their frequencies and compute their significance
  stuffA = dict();
  NA = len(art_wA);
  if NA > 0:
    for com in freqA:
      cs = comm_size[com]
      stuffA[com] = dict()
      L = freqA[com].items()
      L.sort(cmpval)
      for i in range(min(10,len(L))):
        auth = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaA[auth] * 1.0 / NA 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffA[com][i] = [auth.replace('[anonymous]','anonymous').replace('[no Author name AVAILABLE]','anonymous'), f*100, sigma]

  #######
  # INSTITUTIONS
  if verbose: print "....most frequent institutions"
  art_wI = dict(); ## lists article with institutions
  probaI = dict(); ## records the freq of each institution in the whole database
  freqI  = dict(); ## "" in each community of size > thr
  ## the following dict() are used to ensure that we count each pair "article-institution" only once
  probaI_aux = dict(); ## records the freq of each institutions in the whole database
  freqI_aux  = dict(); ## "" in each community of size > thr
  for com in list_comm:  
    freqI[com] = dict()
    freqI_aux[com] = dict()

  # read data
  pl = Utils.Institution()
  pl.read_file(src7)  
  for l in pl.institutions:
      # list the articles with institutions
      if l.id not in art_wI: art_wI[l.id] = ''
      # record the number of occurrence of a given I
      if l.institution not in probaI: 
        probaI[l.institution] = 0; probaI_aux[l.institution] = []
      if l.id not in probaI_aux[l.institution]:
        probaI_aux[l.institution].append( l.id )
        probaI[l.institution] += 1
      # record the occurrence of a given I within each community
      if l.id in list_id: 
        com = partition[l.id]
        if l.institution not in freqI[com]: 
          freqI[com][l.institution] = 0; freqI_aux[com][l.institution] = [] 
        if l.id not in freqI_aux[com][l.institution]:
          freqI_aux[com][l.institution].append( l.id )
          freqI[com][l.institution] += 1
  #
  # extract the 20 most frequent institutions within each community, normalize their frequencies and compute their significance
  stuffI = dict();
  NI = len(art_wI);
  if NI > 0:
    for com in freqI:
      cs = comm_size[com]
      stuffI[com] = dict()
      L = freqI[com].items()
      L.sort(cmpval)
      for i in range(min(20,len(L))):
        inst = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaI[inst] * 1.0 / NI 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffI[com][i] = [inst.replace('&','\\&'), f*100, sigma]

  #######
  # COUNTRIES
  if verbose: print "....most frequent countries"
  art_wC = dict(); ## lists article with countries
  probaC = dict(); ## records the freq of each country in the whole database
  freqC  = dict(); ## "" in each community of size > thr
  ## the following dict() are used to ensure that we count each pair "article-country" only once
  probaC_aux = dict(); ## records the freq of each institutions in the whole database
  freqC_aux  = dict(); ## "" in each community of size > thr
  for com in list_comm:  
    freqC[com] = dict() 
    freqC_aux[com] = dict();

  # read data
  pl = Utils.Country()
  pl.read_file(src6)  
  for l in pl.countries:
      # list the articles with countries
      if l.id not in art_wC: art_wC[l.id] = ''
      # record the number of occurrence of a given C
      if l.country not in probaC: 
        probaC[l.country] = 0; probaC_aux[l.country] = []; 
      if l.id not in probaC_aux[l.country]:
        probaC_aux[l.country].append( l.id ); 
        probaC[l.country] += 1
      # record the occurrence of a given C within each community
      if l.id in list_id: 
        com = partition[l.id]
        if l.country not in freqC[com]: 
          freqC[com][l.country] = 0; freqC_aux[com][l.country] = [];
        if l.id not in freqC_aux[com][l.country]:
          freqC_aux[com][l.country].append( l.id ); 
          freqC[com][l.country] += 1

  # extract the 10 most frequent countries within each community, normalize their frequencies and compute their significance
  stuffC = dict();
  NC = len(art_wC);
  if NC > 0:
    for com in freqC:
      cs = comm_size[com]
      stuffC[com] = dict()
      L = freqC[com].items()
      L.sort(cmpval)
      for i in range(min(10,len(L))):
        coun = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaC[coun] * 1.0 / NC 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffC[com][i] = [coun.replace('&','\\&'), f*100, sigma]

  #######
  # REFERENCES
  if verbose: print "....most frequent references"
  art_wR = dict(); ## lists article with refs
  probaR = dict(); ## records the freq of each ref in the whole database
  freqR  = dict(); ## "" in each community of size > thr
  for com in list_comm: freqR[com] = dict()

  # read data
  pl = Utils.Ref()
  pl.read_file(src5)  
  for l in pl.refs:
      if str(l.volume) !='0':
        reference = l.firstAU + ", " + str(l.year) + ", " + l.journal + " (" + str(l.volume) + "), " + str(l.page)
      else:
        reference = l.firstAU + ", " + str(l.year) + ", " + l.journal
      # list the articles with references
      if l.id not in art_wR: art_wR[l.id] = ''
      # record the number of occurrence of a given R
      if reference not in probaR: probaR[reference] = 0
      probaR[reference] += 1
      # record the occurrence of a given R within each community
      if l.id in list_id: 
        com = partition[l.id]
        if reference not in freqR[com]: freqR[com][reference] = 0
        freqR[com][reference] += 1
  #
  # extract the 20 most frequent refs within each community, normalize their frequencies and compute their significance
  stuffR = dict();
  NR = len(art_wR);
  if NR > 0:
    for com in freqR:
      cs = comm_size[com]
      stuffR[com] = dict()
      L = freqR[com].items()
      L.sort(cmpval)
      for i in range(min(20,len(L))):
        ref = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaR[ref] * 1.0 / NR 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffR[com][i] = [ref.replace('[','').replace(']','').replace('&','\\&'), f*100, sigma]

  #######
  # REFS JOURNALS
  if verbose: print "....most frequent refs journals"
  art_wRJ = dict(); ## lists article with refj
  probaRJ = dict(); ## records the freq of each refj in the whole database
  freqRJ  = dict(); ## "" in each community of size > thr
  ## the following dict() are used to ensure that we count each pair "article-refj" only once
  probaRJ_aux = dict(); ## records the freq of each institutions in the whole database
  freqRJ_aux  = dict(); ## "" in each community of size > thr
  for com in list_comm:  
    freqRJ[com] = dict()
    freqRJ_aux[com] = dict()

  # read data
  pl = Utils.Ref()
  pl.read_file(src5)  
  for l in pl.refs:
      # list the articles with ref journals
      if l.id not in art_wRJ: art_wRJ[l.id] = ''
      # record the number of occurrence of a given RJ
      if l.journal not in probaRJ: 
        probaRJ[l.journal] = 0; probaRJ_aux[l.journal] = dict()
      if l.id not in probaRJ_aux[l.journal]:
        probaRJ_aux[l.journal][l.id] = ''
        probaRJ[l.journal] += 1
      # record the occurrence of a given RJ within each community
      if l.id in list_id: 
        com = partition[l.id]
        if l.journal not in freqRJ[com]: 
          freqRJ[com][l.journal] = 0
          freqRJ_aux[com][l.journal] = dict()
        if l.id not in freqRJ_aux[com][l.journal]:
          freqRJ_aux[com][l.journal][l.id] = '' 
          freqRJ[com][l.journal] += 1
  #
  #
  # extract the 10 most frequent refj within each community, normalize their frequencies and compute their significance
  stuffRJ = dict();
  NRJ = len(art_wRJ);
  if NRJ > 0:
    for com in freqRJ:
      cs = comm_size[com]
      stuffRJ[com] = dict()
      L = freqRJ[com].items()
      L.sort(cmpval)
      for i in range(min(10,len(L))):
        refj = L[i][0]
        f = L[i][1] * 1.0 / cs 
        p = probaRJ[refj] * 1.0 / NRJ 
        if p < 1: sigma = math.sqrt(cs) * (f - p) * 1.0 / math.sqrt(p*(1-p)) 
        else: sigma = 0
        stuffRJ[com][i] = [refj.replace('&','\\&'), f*100, sigma]

  ## END
  #stuffI = replace_accent(stuffI)
  #stuffJ = replace_accent(stuffJ)
  return (stuffK, stuffTK, stuffS, stuffJ, stuffA, stuffI, stuffC, stuffR, stuffRJ)
## ##################################################

def replace_accent(XX):
  for com in XX:
    for i in XX[com]:
      a=1
  return XX

## ##################################################
def direction_and_glue(G,in_dir,partition,thr,verbose,txt):

  #. initialize
  list_nodes = dict();
  comm_size = dict();
  for com in set(partition.values()) :
    list_nodes[com] = [nodes for nodes in partition.keys() if partition[nodes] == com]
    comm_size[com] = len(list_nodes[com])
  #    
  list_comm = [] ## list of communities of size > thr
  list_id = dict(); ## list of id -- articles within the BC network and within a community of size > thr
  for com in comm_size: 
    if comm_size[com] > thr: list_comm.append( com )
  for id_art in partition:
    com = partition[id_art]
    if comm_size[com] > thr: list_id[id_art] = ''   

  #. cluster ref
  if verbose: print"....prep clusters refs table"
  freqR  = dict(); ## records the freq of each ref in each community of size > thr
  norm = dict(); ## tot of previous doc for normalization
  for com in list_comm: 
    freqR[com] = dict()
    norm[com] = 0
  pl = Utils.Ref()
  src  = os.path.join(in_dir, "references.dat")
  pl.read_file(src)  
  for l in pl.refs:
      if str(l.volume) !='0':
        reference = l.firstAU + ", " + str(l.year) + ", " + l.journal + " (" + str(l.volume) + "), " + str(l.page)
      else:
        reference = l.firstAU + ", " + str(l.year) + ", " + l.journal
      # record the occurrence of a given R within each community
      if l.id in list_id: 
        com = partition[l.id]
        if reference not in freqR[com]: freqR[com][reference] = 0
        freqR[com][reference] += 1
  #      
  for com in list_comm:
    for ref in freqR[com]:
      norm[com] += freqR[com][ref]*freqR[com][ref]
    norm[com] = math.sqrt(norm[com]) 
  #
  name = "refs_%s.csv" % (txt)
  dst = os.path.join(in_dir, name)      
  f_r = open(dst,'w')
  f_r.write('com\tref\tn\n')      
  for com in list_comm:
    aux=freqR[com].items()
    aux.sort(cmpval) 
    for kpt in range(min(len(aux),50)):
      f_r.write('%d\t%s\t%d\n' % (com,aux[kpt][0],aux[kpt][1]))      
  f_r.close() 

  #. direction of links btw clusters
  if verbose: print "....prep coupling refs table"
  name = "couplingrefs_%s.csv" % (txt)
  dst = os.path.join(in_dir, name)
  f_cr = open(dst,'w')
  f_cr.write('com1\tcom2\tref\tn12\tn1\tn2\n')
  direction = dict()
  for com1 in list_nodes:
    for com2 in list_nodes:
      size1 = len(list_nodes[com1]); size2 = len(list_nodes[com2]);
      if size1 > thr and size2 > thr and com1 > com2:
        normsh=0; cos1=0; cos2=0;
        shared_ref=[ref for ref in freqR[com2] if ref in freqR[com1]]
        aux=dict()
        for ref in shared_ref:
          cos1 += freqR[com1][ref]*freqR[com1][ref]*freqR[com2][ref]
          cos2 += freqR[com2][ref]*freqR[com1][ref]*freqR[com2][ref]
          normsh += freqR[com1][ref]*freqR[com2][ref]*freqR[com1][ref]*freqR[com2][ref]
          aux[ref] = freqR[com1][ref]*freqR[com2][ref]
        normsh=math.sqrt(normsh) 
        if normsh >0:
          cos1 *= 1.0/(norm[com1]*normsh)
          cos2 *= 1.0/(norm[com2]*normsh)
        direction[(com1,com2)]=[cos1,cos2]
        #.. top 50 coupling ref
        coupling=aux.items()
        coupling.sort(cmpval)
        if cos1>cos2:
          for kpt in range(min(len(coupling),50)):
            ref=coupling[kpt][0]
            f_cr.write('%d\t%d\t%s\t%d\t%d\t%d\n' % (com1, com2, ref, aux[ref], freqR[com1][ref], freqR[com2][ref]))
        else:
          for kpt in range(min(len(coupling),50)):
            ref=coupling[kpt][0]
            f_cr.write('%d\t%d\t%s\t%d\t%d\t%d\n' % (com2, com1, ref, aux[ref], freqR[com2][ref], freqR[com1][ref]))    
  f_cr.close()      

  #. output links with local glue > 0
  list_ref=dict()
  # list_ref['Gardner H, 1983, FRAMES MIND THEORY M']=''
  # list_ref['Kuhn T, 1970, STRUCTURE SCI REVOLU']=''
  # list_ref['Vygotsky L, 1978, MIND SOC DEV HIGHER']=''
  # list_ref['Bandura A, 1986, SOCIAL FDN THOUGHT A']=''
  # list_ref['National council of teachers of MATHEMATICS, 1989, CURR EV STAND SCH MA']=''
  # list_ref['Cohen J, 1983, APPL MULTIPLE REGRES']=''
  # list_ref['Cohen J, 1988, STAT POWER ANAL BEHA']=''
  # list_ref['Reif F, 1965, FUNDAMENTALS STAT TH']=''
  if len(list_ref) >0:
    if verbose: print"....compute local glue"
    for r in list_ref: print r
    name = "links_dir_glue_%s.txt" % (txt)
    dst = os.path.join(in_dir, name)
    f_gephi = open(dst,'w')
    f_gephi.write("edgedef>node1 VARCHAR,node2 VARCHAR,num_links DOUBLE,weight DOUBLE,logweight DOUBLE,cos1 DOUBLE,cos2 DOUBLE")
    for kk in range(len(list_ref)):
      f_gephi.write(',locglue%s DOUBLE' % kk)
    f_gephi.write('\n')  
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
              f_gephi.write("%d,%d,%d,%1.9f,%1.2f,%.5f,%1.5f" % (com1, com2, N, W, 6 + math.log(W)/math.log(10),cos1,cos2))
            else:
              f_gephi.write("%d,%d,%d,%1.9f,%1.2f,%.5f,%1.5f" % (com2, com1, N, W, 6 + math.log(W)/math.log(10),cos2,cos1))  
            #
            for r in list_ref:
              z=0
              if r in freqR[com1] and r in freqR[com2]: z=1 
              f_gephi.write(',%d' % z)
            f_gephi.write('\n')  
    ## ... end
    f_gephi.close() 

  #. compute global glue
  if verbose: print"....compute global glue"
  glue = dict()
  for com1 in list_nodes:
    for com2 in list_nodes:
      size1 = len(list_nodes[com1]); size2 = len(list_nodes[com2]);
      if size1 > thr and size2 > thr and com1 > com2:
        #.. num_links
        N = 0;
        for id1 in list_nodes[com1]:
          for id2 in list_nodes[com2]:
            if id2 in G.edge[id1]: 
              N += 1
        #..add all ref in glue
        if N>0:
          for ref in freqR[com1]:
            if ref in freqR[com2]:
              if ref not in glue: glue[ref]=0
              glue[ref] += freqR[com1][ref]*freqR[com2][ref]*1.0/N
  #. normalize factor
  ZZ=sum(glue.values())
  #. export the 1000 more gluing refs   
  stuff = dict()
  L = glue.items()
  L.sort(cmpval)
  for i in range(min(1000,len(L))):
    r=L[i][0]
    g=L[i][1]*100.0/ZZ
    stuff[i]=[r,g]

  ## END
  return [stuff, direction]

## ##################################################
## ##################################################

def cmpval(x,y):
    if x[1]>y[1]:
        return -1
    elif x[1]==y[1]:
        return 0
    else:
        return 1

## ##################################################
## ##################################################


def main():
# usage: BC.py [-h] [--version] -i DIR -o DIR -p FILE [-v]
#
# optional arguments:
#   -i DIR, --input_dir DIR input directory name
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
          help="output directory name",
          metavar='DIR')
          
  parser.add_argument("-p", "--partition",
          action = "store", dest="partition")

  parser.add_argument("-thr", "--thr", nargs=1, type=int,
          action = "store", dest="thr", default=[10],
          help="",
          metavar='thr')

  parser.add_argument("-v", "--verbose",
          action = "store_true", dest="verbose",
          default = False,
          help="verbose mode [default %(default)s]")

  #Analysis of input parameters
  args = parser.parse_args()
  
  if (not os.path.exists(args.in_dir[0])):
      print "Error: Input directory does not exist: ", args.in_dir[0]
      exit()

  if (not os.path.exists(args.partition[0])):
      print "Error: Input partition does not exist: ", args.partition[0]
      exit()

  if (not os.path.exists(args.out_dir[0])):
      print "Error: Output directory does not exist: ", args.out_dir[0]
      exit()

  ##

  partition = dict();      

  f_in = open(args.partition,'r')
  for line in f_in.readlines():
    foo = line.stript().split('\t')
    if len(foo) == 2:
      partition[int(foo[0])] = int(foo[1])
  f_in.close()
  
  ##
  comm_tables(args.in_dir[0],partition,thr,args.verbose)

  return


    
## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################
