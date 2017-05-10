import sys, os
from statistics import median
import matplotlib.pyplot as plt
from math import sqrt, log
import pandas as pd
from bisect import bisect


if len(sys.argv) < 2:
	print("Give the image filename")
	sys.exit(2);
k = 13 if len(sys.argv) < 3 else int(sys.argv[2]) # number of neighbors taken into account
verbose = len(sys.argv) > 3 and sys.argv[3] == "v"
maxFeatures = 20

stream = os.popen("./src/invariants "+ sys.argv[1]).read()
img = [float(j) for j in stream.split(",")]
maxFeatures = min(maxFeatures, len(img))

data = pd.read_csv("src/dataAll.csv")
N = data.shape[0]

def featName(i):
	return data.columns[i+1]

##################### normalization of features #####################
def normalize(x):
	return float(x) / N

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

################################ kNN ################################
f = open("classes.csv","r")
classes = [c.split(",")[0] for c in f.read().split("\n")]
f.close()

vbOutput = ""

proba = {}
probaMin = 0.1 # equally distributed proba 
probaNeighbor = 0.7 # proba given to k neighbors
probaProximity = 1 - probaMin - probaNeighbor # proba given wrt proximity
unitMin = probaMin / len(classes)
unitNeighbor = probaNeighbor / k
unitProximity = probaProximity / (N*(N+1)/2)

distances.sort()
for i in range(k):
	index = distances[i][1]
	c = data.ix[index]['category']
	vbOutput += "{}*{}\t\t{}\n".format(c, index%15, distances[i][0])
	proba[c] = proba.get(c, unitMin) + unitNeighbor + unitProximity * (N - i)
for i in range(k,N):
	c = data.ix[distances[i][1]]['category']
	proba[c] = proba.get(c, unitMin) + unitProximity * (N - i)

############################## output ###############################
if not verbose:
	txt = "\n".join([str(proba.get(c, unitMin)) for c in classes])
	print(txt)
else:
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
