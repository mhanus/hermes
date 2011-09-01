nx1 = 10
ny1 = 10
nx2 = 20
ny2 = 20
dx1 = 65./nx1;
dy1 = 60./ny1;
dx2 = (133.-65.)/nx2;
dy2 = (140.-60.)/ny2;

s1 = ''
s2 = ''
n1 = 1
n2 = 1

y = 0.
for j in range(0,ny1+ny2):
  if y < 60:
    dy = dy1
  else:
    dy = dy2

  x = 0.
  for i in range(0,nx1+nx2):
    if x < 65:
      dx = dx1
    else:
      dx = dx2
      
    if (i == nx1):
      s1 = s1 + repr(n1) + ' '
      n1 = n1 + 1
      print repr(x+dx2/2.) + ", " + repr(y+dy/2.)
    else:
      s1 = s1 + '0 '
      
    if (j == 0):
      s2 = s2 + repr(n2) + ' '
      n2 = n2+1
      print repr(x+dx/2.) + ", " + repr(y+dy1/2.)
    else:
      s2 = s2 + '0 '
    
    x = x+dx      
  
  y = y+dy
  
  s1 = s1 + '\n'
  s2 = s2 + '\n'
    
print s1
print "________________________________________________________________________"
print
print s2
