#! /usr/bin/env python

""" 
   Author : Sebastian Grauwin (http://www.sebastian-grauwin.com/)
   Copyright (C) 2012
   All rights reserved.
   BSD license.
   .... If you are using these scripts, please cite our "Scientometrics" paper:
   .... S Grauwin, P Jensen, Mapping Scientific Institutions. Scientometrics 89(3), 943-954 (2011)
"""

# usage: parser.py -i DIR -o DIR [-v]
# 

import os
import sys
import glob
import re
import numpy
import math
import time
import argparse
#import Utils.Utils as Utils
import UtilsBeta.Utils as Utils
from fnmatch import fnmatch

## ##################################################
# a dictionnary of common_words used when extracting title words
common_words = ['a', 'able', 'about', 'across', 'after', 'all', 'almost', 'also', 'am', 'among', 'an', 'and', 'any', 'are', 'as', 'at', 'be', 'because', 'been', 'but', 'by', 'can', 'cannot', 'could', 'dear', 'did', 'do', 'does', 'either', 'else', 'ever', 'every', 'for', 'from', 'get', 'got', 'had', 'has', 'have', 'he', 'her', 'hers', 'him', 'his', 'how', 'however', 'i', 'if', 'in', 'into', 'is', 'it', 'its', 'just', 'least', 'let', 'like', 'likely', 'may', 'me', 'might', 'most', 'must', 'my', 'neither', 'no', 'nor', 'not', 'of', 'off', 'often', 'on', 'only', 'or', 'other', 'our', 'own', 'rather', 'said', 'say', 'says', 'she', 'should', 'since', 'so', 'some', 'than', 'that', 'the', 'their', 'them', 'then', 'there', 'these', 'they', 'this', 'tis', 'to', 'too', 'twas', 'us', 'wants', 'was', 'we', 'were', 'what', 'when', 'where', 'which', 'while', 'who', 'whom', 'why', 'will', 'with', 'would', 'yet', 'you', 'your']
french_words = ['un','une','au','aux','entre','a','le','la','les','du','de','des','mais','ou','et','dans','avec','en','sur','sous','avant','apres','vers','par','pendant','depuis','pour','chez','est','ont']
common_words = common_words + french_words;
punctuation = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '.', '/', ':', ';', '<', '=', '>', '?', '@', '[', '\\', ']', '^', '_', '`', '{', '|', '}', '~', ' - ']
# I kept the '-' out of this list


## ##################################################
## ##################################################

def biblio_parser(in_dir,out_dir,database,verbose):

  ## INITIALIZATION
  t1=time.time()
  #.. detecting raw files
  if database == 'wos': pattern = "*.txt"
  if database == 'scopus': pattern = "*.csv"
  if verbose: print "..Analysing files %s/%s" % (in_dir,pattern)
  srclst=[]
  for path, subdirs, files in os.walk(in_dir):
    for name in files:
        if fnmatch(name, pattern):
            srclst.append( os.path.join(path, name))
  print "....%d '%s' files detected" % (len(srclst),pattern)
  #.. filtering according to publication year?
  myfilter=raw_input("....do you want to put some filter on the publication year? (y/n): ")
  if (myfilter=='y'):
    confirm='n';
    while (confirm!='y'):
      y_min=input("..Publication year must be >=")
      y_max=input("..Publication year must be <=")
      confirm=raw_input("..Keep only papers with publication year between %d and %d. Confirm? (y/n): " % (y_min,y_max))
  else:
    y_min=1500
    y_max=3000
  #.. prep empty parsed files
  dst1  = os.path.join(out_dir, "articles.dat") 
  f_articles = open(dst1,'w')
  dst2  = os.path.join(out_dir, "authors.dat")
  f_authors = open(dst2,'w')
  dst3  = os.path.join(out_dir, "keywords.dat")
  f_keywords = open(dst3,'w')
  dst4  = os.path.join(out_dir, "subjects.dat")
  f_subjects = open(dst4,'w')
  dst5  = os.path.join(out_dir, "references.dat")
  f_refs = open(dst5,'w')
  dst6  = os.path.join(out_dir, "countries.dat")
  f_countries = open(dst6,'w')
  dst6b  = os.path.join(out_dir, "cities.dat")
  f_cities = open(dst6b,'w')
  dst7  = os.path.join(out_dir, "institutions.dat")
  f_institutions = open(dst7,'w')
  dst8  = os.path.join(out_dir, "abstracts.dat")
  f_abstracts = open(dst8,'w')  
  #.. some parameters to count errors 
  kompt_refs = 0
  kompt_corrupt_refs = 0
  #.. "id" will be the new id for articles in the dataset, "unique_ids" tracks the WOS or SCOPUS unique id that we use to remove duplicates
  id = int(-1)
  UNIQUE_IDS = dict() 

  ## SCOPUS JOURNAL CATEGORIES
  #... scopus export data do not contain any subject category. Here we upload an official scopus list from february 2015 listing the SUBJCAT each journal correspoind to (multiple category per journal possible)
  if (database =='scopus'):
    if (verbose): print "..upload Scopus journal categories from auxiliary file"
    keep = raw_input("....keep Scopus Categories (Ncat=26) or Sub-categories (Nsubcat=308)? (C/S):")
    journalCATS=dict()
    code_cat=dict()
    script_dir = os.path.dirname(__file__)
    #.. read cat_catcode file
    filename = os.path.join(script_dir, 'Utils/scopus_subcat_codes.txt')
    with open(filename,'r') as fd:
      for line in fd.readlines():
        line = line.strip('\n')
        foo=line.split('\t')
        code_cat[foo[1]]=foo[0]
    #.. read journal_catcode file
    filename = os.path.join(script_dir, 'Utils/scopus_journals_subcat.txt')
    with open(filename,'r') as fd:
      for line in fd.readlines():
        line = line.strip('\n')
        foo=line.split('\t')
        jrn=foo[0]
        bar=foo[1].replace('; ',';').split(';')
        journalCATS[jrn]=[]
        for i in range(len(bar)):
          if len(bar[i])>0:
            if keep=='S': 
              cat=code_cat[bar[i]]
            if keep=='C':  
              cat=code_cat[bar[i][0:2]+'00'].replace('(all)','')
            journalCATS[jrn].append(cat)  
    for jrn in journalCATS: journalCATS[jrn]=list(set(journalCATS[jrn])) 
    del code_cat;     
    # create dict of journals / sources not in that category file        
    journalNOT=dict()

  ## TREAT DATA
  # will parse info about each article and write clean version in output files
  for src in srclst:
      pl = Utils.ArticleList()
      pl.read_file(src,database)
      if verbose: print "..processing %d articles in file %s" % (len(pl.articles), src)  
      if (len(pl.articles) > 0):
          for article in pl.articles:
            if article.UT not in UNIQUE_IDS and (article.PY >= y_min and article.PY <= y_max):
              UNIQUE_IDS[article.UT] = ''
              id = id + 1
              #article 
              foo = article.AU.split('; ')
              firstAU = foo[0].replace(',','')
              if (article.J9==''): article.J9='[]';
              if (article.DT==''): article.DT='[]';
              f_articles.write("%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (id,firstAU,article.PY,article.J9,article.VL,article.BP,article.DI,article.PT,article.DT,article.LA,article.TC,article.TI,article.UT))
              #authors
              if(article.AU != ""): 
                  foo = article.AU.split('; ')
                  for i in range(len(foo)):
                      #... check the upper or lower format of letter to have the author name in format "Grauwin S"
                      foo[i] = foo[i].replace(',','')
                      aux1 = foo[i].rfind(' ')
                      aux2 = len(foo[i])
                      foobar = foo[i].lower().capitalize()
                      if aux1 > 0: 
                          s1 = foobar[aux1:aux2]
                          s2 = s1.upper() 
                          foobar = foobar.replace(s1,s2)
                      aux = foobar.find('-')
                      if aux > 0: 
                          bar1 = foobar[aux:aux+2]
                          bar2 = '-' + foobar[aux+1].upper()
                          foobar = foobar.replace(bar1,bar2)
                      aux = foobar.find(' ')
                      if aux > 0: 
                          bar1 = foobar[aux:aux+2]
                          bar2 = ' ' + foobar[aux+1].upper()
                          foobar = foobar.replace(bar1,bar2)
                      f_authors.write("%d\t%d\t%s\n" % (id,i,foobar))
              #keywords
              if(article.DE != ""):
                  #.. author keywords
                  foo = article.DE.split('; ')
                  for f in foo:
                      f_keywords.write("%d\tAK\t%s\n" % (id,f.upper()))
              if(article.ID != ""):
                  #.. WOS or SCOPUS keywords 
                  foo = article.ID.split('; ')
                  for f in foo:
                      f=re.sub('[\(\)]','', f) # (remove stuff in parentheses - happen sometimes in scopus)
                      f_keywords.write("%d\tIK\t%s\n" % (id,f.upper()))
              if(article.TI != ""):
                  #.. title words (excluding the common words listed on the top of this file)
                  foo = article.TI
                  #... remove ponctuations !"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~
                  for p in punctuation: foo = foo.replace(p,'')
                  foo = foo.split(' ')
                  for f in foo:
                    bar = f.lower()
                    if bar not in common_words and len(bar)>0:
                      f_keywords.write("%d\tTK\t%s\n" % (id, bar.upper()))
              #subjects
              if database=='wos':
                if(article.WC != ""):
                  foo = article.WC.split('; ')
                  for cat in foo: f_subjects.write("%d\t%s\n" % (id,cat))
              if database=='scopus':
                if article.SO in journalCATS:
                  foo=journalCATS[article.SO]
                  for cat in foo: f_subjects.write("%d\t%s\n" % (id,cat))
                else: 
                  if article.SO not in journalNOT: journalNOT[article.SO]=0
                  journalNOT[article.SO]+=1;  
              # #references
              if(article.CR != ""):
                   foo = article.CR.split('; ')
                   for i in range(len(foo)):
                       ref=Utils.Ref()
                       ref.parse_ref(foo[i],database)
                       kompt_refs += 1 
                       if(ref.year > 0): 
                           f_refs.write("%d\t%s\t%d\t%s\t%s\t%s\n" % (id,ref.firstAU,ref.year,ref.journal,ref.volume,ref.page))
                       if(ref.year == 0):
                         #print foo[i]
                         kompt_corrupt_refs += 1  
              #countries / cities / institutions
              if(article.C1 != ""):
                  adresse = article.C1
                  aux1 = adresse.find('[')
                  aux2 = adresse.find(']')
                  while (aux1 < aux2) and aux1 > -1:
                      aux = adresse[aux1:min(aux2+2,len(adresse))]
                      adresse = adresse.replace(aux,'')
                      aux1 = adresse.find('[')
                      aux2 = adresse.find(']')
                  foo = adresse.split('; ')
                  for i in range(len(foo)):
                      #... "institutions" will keep everything listed between commas in the adresses, except from cities and countries
                      foo[i] = foo[i].replace(', ',',')
                      bar = foo[i].split(',') 
                      ll = len(bar)
                      for j in range(ll - 2):
                        f_institutions.write("%d\t%d\t%s\n" % (id,i,bar[j]))
                      #... country
                      country = bar[ll-1].replace('.','').replace(';','').upper()
                      lll=len(bar[ll-1])
                      #....  if pb, indicate country X
                      if lll<2: country='X'
                      #.... put all USA states together under the label 'USA'
                      usa_states=['AL','AK','AZ','AR','CA','NC','SC','CO','CT','ND','SD','DE','FL','GA','HI','ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS','MO','MT','NE','NV','NH','NJ','NM','NY','OH','OK','PA','RI','TN','TX','UT','VT','VA','WV','WA','WI','WY','DC'];
                      usa_states2=[f+' ' for f in usa_states];
                      if (country[lll-3:lll] == 'USA' or country in usa_states or country[0:3] in usa_states2): country = 'USA'
                      #.... put England, Wales, North Ireland, Scotland in UK
                      if country in ['ENGLAND','WALES','NORTH IRELAND','SCOTLAND','UNITED KINGDOM']: country='UK'
                      if country not in ['USA','UK']: 
                        country=" ".join(w.capitalize() for w in country.lower().split())
                      if (database =='scopus' and country == 'USA'):country='United States'
                      if (database =='scopus' and country == 'UK'):country='United Kingdom'
                      f_countries.write("%d\t%d\t%s\n" % (id,i,country))
                      #... cities
                      if country == 'France':
                        city=bar[ll-2];
                        for nb in range(10): 
                          city=city.replace(str(nb),'');
                        city=city.replace('F-','');
                        city=city.replace('FR-','');
                        city=city.replace('FR -','');
                        city=city.lower().replace('cedex','');
                        city=" ".join(w.capitalize() for w in city.split(' ') if len(w)>0)
                        #if len(city) >0:
                        # if city[-1]==' ': city=city[0:-1];
                        # if city[0]==' ': city=city[1:];
                        city=city.lower().capitalize()  
                        if len(city) >0:
                          f_cities.write("%d\t%d\t%s\n" % (id,i,city))
              # abstract (the 2 conditions refer to wos and scopus way to indicate the absence of abstract)
              if(article.AB != "") and (article.AB != "[No abstract available]"):
                f_abstracts.write("%d\t%s\n" % (id,article.AB))
      #.. delete the stuff dealing from the raw file just treated from memory        
      del pl

  ## END
  #... how many papers?
  print("..%d parsed articles in total") % (id + 1) 
  #... error in refs?
  if kompt_refs > 0: print("..%d inadequate refs out of %d (%f%%) have been rejected by this parsing process (no publication year, unpublished, ...) ") % (kompt_corrupt_refs, kompt_refs, (100.0 * kompt_corrupt_refs) / kompt_refs)
  else: print '..no references found in your dataset. Check whether you extracted your data properly!'
  #... journal with unknown cat?
  if database =='scopus':
    print "..%d publication sources with unknown subject category" % (len(journalNOT))
    if len(journalNOT)>0:
      foo=len([jr for jr in journalNOT if journalNOT[jr]>=5])
      print "..%d of them with more than 5 papers in the corpus" % (foo)
      filename = 'scopus_sources_without_categories.dat'
      print "...all these publication sources are listed in %s" % filename
      with open(filename,'r') as out:
        out.write('Publication source\tNpapers\n')
        fb=journalNOT.items()
        fb.sort(cmpval)
        for elt in fb:
          out.write("%s\t%d\n" % (elt[0],elt[1]))
  #... close output files
  f_articles.close()
  f_authors.close()
  f_keywords.close()
  f_subjects.close()
  f_refs.close()
  f_countries.close()
  f_cities.close()
  f_institutions.close()
  f_abstracts.close()

  t2=time.time()
  print 'total time needed: %ds' % (t2-t1)
  return

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
## ##################################################

def main():
# usage: parser.py [-h] [--version] -i DIR -o DIR [-v]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   --version             show program's version number and exit
#   -o DIR, --output_dir DIR
#                         output directory name
#   -i DIR, --input_dir DIR
#                         input directory name
  # Parse line options.
  # Try to have always the same input options
  parser = argparse.ArgumentParser(description = 'parser')

  parser.add_argument('--version', action='version', version='%(prog)s 1.1')
  
  parser.add_argument("-i", "--input_dir", nargs=1, required=True,
          action = "store", dest="in_dir",
          help="input directory name",
          metavar='DIR')
          
  parser.add_argument("-o", "--output_dir", nargs=1, required=True,
          action = "store", dest="out_dir",
          help="output directory name",
          metavar='DIR')

  parser.add_argument("-d", "--database",
          action = "store", dest="database",
          default = 'wos',
          help="database [default %(default)s]",
          metavar='string')  

  parser.add_argument("-v", "--verbose",
          action = "store_true", dest="verbose",
          default = False,
          help="verbose mode [default %(default)s]")

  #Analysis of input parameters
  args = parser.parse_args()
  
  if (not os.path.exists(args.in_dir[0])):
    print "Error: Input directory does not exist: ", args.in_dir[0]
    exit()

  if (not os.path.exists(args.out_dir[0])):
    print "Error: Output directory does not exist: ", args.out_dir[0]
    exit()

  if args.database not in ['wos','scopus']:
    print "Error: database must be either 'wos' or 'scopus'"
    exit()

  ##      

  biblio_parser(args.in_dir[0],args.out_dir[0],args.database,args.verbose)

  return


    
## ##################################################
## ##################################################
## ##################################################

if __name__ == "__main__":
    main()

## ##################################################
## ##################################################
## ##################################################

