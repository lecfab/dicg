import sys, os
from statistics import median
import matplotlib.pyplot as plt
from math import sqrt, log


if len(sys.argv) < 2:
	print("Give the image filename")
	sys.exit(2);
k = 11 if len(sys.argv) < 3 else int(sys.argv[2]) # number of neighbors taken into account
verbose = len(sys.argv) > 3 and sys.argv[3] == "v"
maxFeatures = 2

stream = os.popen("./src/invariants "+ sys.argv[1]).read()
img = [0] + [float(j) for j in stream.split(",")]
maxFeatures = min(maxFeatures, len(img)-1)

data = open("src/data_noScale.csv","r")
lines = data.read()
data.close()
lines = list(filter(lambda x: x != '', lines.split("\n")))[1:]

##################### normalization of features #####################
def apply(x):
	return float(x)

imgX = apply(img[1])
imgY = apply(img[2])

featSum = [0] * 10
for line in lines:
	v = line.split(",")
	for i in range(maxFeatures):
		featSum[i+1] += apply(v[i+1])
featAvg = [x / len(lines) for x in featSum]

####################### distances computation #######################
distances = []
metrics = [1] * 10
for line in lines:
	v = line.split(",")
	d = 0
	x = []
	for i in range(1, maxFeatures+1):
		x.append(float(v[i]))
		d += metrics[i] * ((apply(v[i]) - apply(img[i])) / featAvg[i]) ** 2
	distances.append((d,v[0], x[0], x[1]))
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
	c = distances[i][1]
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
		h = hash(c) % (len(colors)*len(symbols))
		hS = int(h / len(colors)) # choose symbol
		hC = h % len(colors) # choose color
		plt.plot(x, (y), colors[hC] + symbols[hS])
		#plt.annotate(c, (x, (y)))
	plt.show()
else:
	print(txt)