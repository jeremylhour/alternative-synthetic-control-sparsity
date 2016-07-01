### Citations of Abadie et al (2010)
### Jeremy L Hour
### 30 juin 2016

import numpy as np

f = open("E:/BEAST/biblioanalysis/abadie2010.txt","r")
d = []
for line in f:
  line = line.strip()
  columns = line.split()
  d.append(int(columns[1]))
f.close()
d=np.array(d)
