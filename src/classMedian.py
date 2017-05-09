from math import log

file = open("inv.csv","r")
lines = file.read()
lines = list(filter(lambda x: x != '', lines.split("\n")))
values = [0]*10
nbFeat = 0
classes = {}

for line in lines:
	v = list(filter(lambda x: x!='', line.split(",")))
	nbFeat = len(v)
	c = v[0]
	if not c in classes:
		classes[c] = [0] * 10
	classes[c][0] += 1
	for i in range(1, nbFeat):
		classes[c][i] += float(v[i])*float(v[i])
file.close()

csv = ""
for c in classes:
	csv += c + ":" + str(classes[c][0])
	for i in range(1, nbFeat):
		classes[c][i] /= classes[c][0]
		csv += ","+"{0:.0f}".format(log(classes[c][i]))
	csv += "\n"

output = open("invClasses.csv", "w")
output.write(csv)
output.close()