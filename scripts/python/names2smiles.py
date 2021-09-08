#names2smiles input output
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
