from math import log

file = open("inv.csv","r")
lines = file.read()
lines = list(filter(lambda x: x != '', lines.split("\n")))
values = [0]*10
total = 0
for line in lines:
	v = list(filter(lambda x: x!='', line.split(",")))
	for i in range(1, len(v)):
		values[i] += float(v[i])*float(v[i])
	total += 1

csv = ""
for line in lines:
	v = list(filter(lambda x: x!='', line.split(",")))
	csv += v[0]
	for i in range(1, len(v)):
		# csv += ","+"{0:f}".format((float(v[i])*float(v[i])*total / values[i]))
		vv = float(v[i])*float(v[i])*total / values[i]
		if(vv != 0): vv = log(vv)
		csv += ","+str(int(vv))
	csv += "\n"
file.close()
output = open("inv2.csv", "w")
output.write(csv)
output.close()