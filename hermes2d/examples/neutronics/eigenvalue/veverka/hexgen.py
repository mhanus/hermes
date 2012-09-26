import numpy, math, pylab

pts_in_layer = [22,11,11,10,10,9,8,8,7,6,6,6,5,5,4,3,3,2,1,1]

nodes = []
h = 14.7
e = h/math.sqrt(3.0)
angle = math.pi/6.0;

for j in range(0,pts_in_layer[0]):
	nodes.append([j*h*0.5,0])

for i in range(1,20):
	y = i*e*math.sin(angle)
	for j in range(0,pts_in_layer[i]):
		x = i*e*math.cos(angle) + j*h
		
		nodes.append([x,y])
		
nodes_array = numpy.array(nodes)

print nodes_array
print nodes_array.shape[0]

numpy.savetxt("core.mesh", nodes_array, fmt='[%1.15g, %1.15g],')


pylab.plot(nodes_array[:,0], nodes_array[:,1],' ')
for i in range(0,nodes_array.shape[0]):
	pylab.text(nodes_array[i,0], nodes_array[i,1], `i`, fontsize=12, horizontalalignment='center', verticalalignment='center')

pylab.axis('equal')
pylab.ylim(ymin=-e)	
pylab.xlim(xmin=-h)

pylab.show()
