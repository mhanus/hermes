import numpy, math, pylab

pts_in_layer = [11,10,10,9,9,8,6,6,5,3,1]

nodes = []
h = 19
angle = math.pi/3.0;

for i in range(0,11):
	y = i*h*math.sin(angle)
	for j in range(0,pts_in_layer[i]):
		x = i*h*math.cos(angle) + j*h
		
		nodes.append([x,y])
		
nodes_array = numpy.array(nodes)

print nodes_array
print nodes_array.shape

numpy.savetxt("core.mesh", nodes_array, fmt='[%1.15g, %1.15g],')

pylab.plot(nodes_array[:,0], nodes_array[:,1],'.')
pylab.show()
