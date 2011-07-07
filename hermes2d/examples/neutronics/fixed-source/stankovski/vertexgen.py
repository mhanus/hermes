import numpy, math, pylab

a = 0.625
r = 0.450
s = [33, 26, 15, 0]

eps = 0.5*math.pi/4.

nodes = []

for i in range(0,8):
	y = i*a
	for j in range(i,8):
		x = j*a
		nodes.append([x,y])

for i in range(4,0,-1):
  id = s[i-1]
  
  if i == 4:
    start = 0.
    end = math.pi/4. + eps
  else:
    start = 5./4.*math.pi
    end = 9./4.*math.pi + eps
      
  for j in range (0,i):
    angle = start
    while angle < end:
      x = nodes[id][0] + r*math.cos(angle)
      y = nodes[id][1] + r*math.sin(angle)
      nodes.append([x,y])
      angle = angle + math.pi/4.
      
    id = id + 2
    
    if i == 4:
      end = math.pi + eps
    else:
      start = 0.
      end = 2*math.pi
    		
nodes_array = numpy.array(nodes)

print nodes_array
print nodes_array.shape

numpy.savetxt("domain.mesh", nodes_array, fmt='[%1.15g, %1.15g],')

pylab.plot(nodes_array[:,0], nodes_array[:,1],'.')
pylab.show()
