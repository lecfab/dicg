import sys, os
from statistics import median
import matplotlib.pyplot as plt
from math import sqrt, log
import pandas as pd
from bisect import bisect


if len(sys.argv) < 2:
	print("Give the image filename")
	sys.exit(2);
k = 90 if len(sys.argv) < 3 else int(sys.argv[2]) # number of neighbors taken into account
verbose = len(sys.argv) > 3 and sys.argv[3] == "v"
maxFeatures = 2

stream = os.popen("./src/invariants "+ sys.argv[1]).read()
img = [float(j) for j in stream.split(",")]
maxFeatures = min(maxFeatures, len(img))

data = pd.read_csv("src/data.csv")

def featName(i):
	return data.columns[i+1]

##################### normalization of features #####################
def normalize(x):
	return float(x) / data.shape[0]

def signature(values):
	s = []
	for i in range(maxFeatures):
		s.append(normalize(values[featName(i)]))
	return s

featA, featB = 0,1 # id's of columns that will be plotted
imgSig = []
for feature in range(maxFeatures):
	dataSorted = data.sort_values(by=featName(feature)).iterrows()
	v = bisect([values[featName(feature)] for i,values in dataSorted], img[feature]) - 0.5
	imgSig.append(normalize(v))

####################### distances computation #######################
def distance(sig1, sig2):
	return sum([abs(sig1[i]-sig2[i])**0.5 for i in range(len(sig1))])

distances = []
for i,values in data.rank().iterrows():
	sig = signature(values)
	d = distance(sig, imgSig)
	distances.append((d, i, sig[featA], sig[featB]))
distances.sort()

################################ kNN ################################
f = open("classes.csv","r")
classes = [c.split(",")[0] for c in f.read().split("\n")]
f.close()

proba = {}
unitary_proba = min(1 / 1000, 0.5 / len(classes)) # minimum proba of a class
neighbors_proba = (1 - len(classes) * unitary_proba) / (k*(k+1)/2) # increase of proba due to neighborhood
vbOutput = ""

for i in range(k):
	index = distances[i][1]
	c = data.ix[index]['category']
	vbOutput += "{}*{}\t\t{}\n".format(c, index%15, distances[i][0])
	# 		   current proba 	+   better proba for closer neighbors
	proba[c] = proba.get(c, unitary_proba) + neighbors_proba * (k - i)

############################## output ###############################
txt = "\n".join([str(proba.get(c, unitary_proba)) for c in classes])
if verbose:
	print(vbOutput)
	colors = "rgbcmyk"
	symbols = "ov^<>1234sp*hH+xDd"
	imgX, imgY = imgSig[featA], imgSig[featB]
	
	plt.axvline(imgX)
	plt.axhline(imgY)
	print([imgX, imgY])
	
	for i in range(k):
		d, c, x, y = distances[i]
		c = data.ix[c]['category']
		h = hash(c) % (len(colors)*len(symbols))
		hS = int(h / len(colors)) # choose symbol
		hC = h % len(colors) # choose color
		plt.plot(x, y, colors[hC] + symbols[hS])
	plt.get_current_fig_manager().window.showMaximized()
	plt.show()
else:
	oouou = open("oouou.txt","w")
	oouou.write(txt)
	print(txt)