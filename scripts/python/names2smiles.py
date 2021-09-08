#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#./name2smiles [input][output]
#Created on Mon Jan 20 12:55:46 2020
#To retireve SMILES based on compounds names or cas number from national cancer institute repositories
import sys
from urllib.request import urlopen
def CIRconvert(ids):
  try:
  url = 'http://cactus.nci.nih.gov/chemical/structure/' + ids + '/smiles'
ans = urlopen(url).read().decode('utf8')
return ans
except:
  return ids + 'Did not work'
import codecs
f_out = codecs.open(sys.argv[2],mode='a',encoding='utf-8')
ids_file=open(sys.argv[1],'r')

for line in ids_file.readlines():
  line=line.replace('\n', '')
print(str(line))
print(line,'\t',CIRconvert(line))
res=line+'\t'+CIRconvert(line)
f_out.write(str(res))
f_out.write('\n')
f_out.close()
