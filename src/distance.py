import sys, os
from statistics import median
import matplotlib.pyplot as plt
from math import sqrt, log
import pandas as pd
from bisect import bisect


if len(sys.argv) < 2:
	print("Give the image filename")
	sys.exit(2);
k = 11 if len(sys.argv) < 3 else int(sys.argv[2]) # number of neighbors taken into account
verbose = len(sys.argv) > 3 and sys.argv[3] == "v"
maxFeatures = 2 # 10

stream = os.popen("./src/invariants "+ sys.argv[1]).read()
img = [0] + [float(j) for j in stream.split(",")]
maxFeatures = min(maxFeatures, len(img)-1)

data = pd.read_csv("src/data.csv")

##################### normalization of features #####################
def apply(x):
	return float(x)

def position(a, b):
	x = float(a)/data.shape[0]
	y = float(b)/data.shape[0]
	y = 0 if a == 0 else float(b)/a
	return x, y

la = data.sort_values(by='a')
lb = data.sort_values(by='b')

imgX = bisect([j['a'] for i, j in la.iterrows()], img[1]) - 0.5
imgY = bisect([j['b'] for i, j in lb.iterrows()], img[2]) - 0.5
imgX, imgY = position(imgX, imgY)

####################### distances computation #######################
distances = []
for i,j in data.rank().iterrows():
	x, y = position(j['a'], j['b'])
	d = (x - imgX) + (y - imgY)
	distances.append((d, i, x, y))
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
	c = data.ix[distances[i][1]]['category']
	vbOutput += "{}\t\t{}\n".format(c, distances[i][0])
	# 		   current proba 	+   better proba for closer neighbors
	proba[c] = proba.get(c, unitary_proba) + neighbors_proba * (k - i)

############################## output ###############################
txt = "\n".join([str(proba.get(c, unitary_proba)) for c in classes])
if verbose:
	print(vbOutput)
	colors = "rgbcmyk"
	symbols = "ov^<>1234sp*hH+xDd"
	plt.axhline(imgY)
	plt.axvline(imgX)
	print([imgX, imgY])
	for i in range(k):
		d, c, x, y = distances[i]
		c = data.ix[c]['category']
		h = hash(c) % (len(colors)*len(symbols))
		hS = int(h / len(colors)) # choose symbol
		hC = h % len(colors) # choose color
		plt.plot(x, y, colors[hC] + symbols[hS])
		#plt.annotate(c, (x, (y)))
	plt.show()
else:
	oouou = open("oouou.txt","w")
	oouou.write(txt)
	print(txt)