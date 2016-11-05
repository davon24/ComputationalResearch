from numpy import*


graphtype = "C:/Users/wdavon24/BarSiteV2/file.txt"
 
data = genfromtxt(graphtype, skiprows = 2, dtype=None, unpack=True) #file name with data points
slope = 0
offset = 0
i = 0
for row in data:
	slope += row[1]/row[0]
	i+=1

slope /= i
print slope

for row in data:
	offset += row[1] - slope*row[0]

offset /= i
print offset