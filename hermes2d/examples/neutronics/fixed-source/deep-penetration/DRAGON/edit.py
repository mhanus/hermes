# Number of mesh subdivisions in each of the two regions.
nx1 = 16
nx2 = 15
nx3 = 15
nx4 = 40
ny1 = 16
ny2 = 10
ny3 = 10
ny4 = 50

ax1 = 50.
ax2 = 65.
ax3 = 80.
ax4 = 133.
ay1 = 50.
ay2 = 60.
ay3 = 70.
ay4 = 140. 

# Length of subdivided mesh intervals.
dx1 = ax1 / nx1
dx2 = (ax2 - ax1) / nx2
dx3 = (ax3 - ax2) / nx3
dx4 = (ax4 - ax3) / nx4
dy1 = ay1 / ny1
dy2 = (ay2 - ay1) / ny2
dy3 = (ay3 - ay2) / ny3
dy4 = (ay4 - ay3) / ny4

i1 = nx1 + nx2

# Resulting strings with edit regions for the two tallies.
s1 = ''
s2 = ''
# Gradually increased numbers representing merge regions in the tallies.
n1 = 1
n2 = 1
# Column counters for wrapping after 68th column.
col1 = 1
col2 = 1

y = 0.
for j in range(0,ny1+ny2+ny3+ny4):
  if y < ay1:
    dy = dy1
  elif y >= ay1 and y < ay2:
    dy = dy2
  elif y >= ay2 and y < ay3:
    dy = dy3
  elif y >= ay3 and y < ay4:
    dy = dy4

  x = 0.
  for i in range(0,nx1+nx2+nx3+nx4):
    if x < ax1:
      dx = dx1
    elif x >= ax1 and x < ax2:
      dx = dx2
    elif x >= ax2 and x < ax3:
      dx = dx3
    elif x >= ax3 and x < ax4:
      dx = dx4  
    
    # First tally (flux at x=65cm).  
    if (i == i1):
      s1 = s1 + str(n1)
      n1 = n1 + 1
      print str(x+dx3/2.) + ", " + str(y+dy/2.)
    else:
      s1 = s1 + '0'
      
    # Second tally (flux at y=0cm).
    if (j == 0):
      s2 = s2 + str(n2)
      n2 = n2 + 1
      print str(x+dx/2.) + ", " + str(y+dy1/2.)
    else:
      s2 = s2 + '0'
     
    # Append a space and either the merge-region index or zero (indicating that
    # the current mesh region is to be ignored in this tally).
    s1 = s1 + ' '
    col1 = col1 + 1 + (len(str(n1-1)) if i == i1 else 1)
    
    # Should the next entry overflow the 68th column, remove the trailing space
    # and break the line.
    if (col1 + (len(str(n1)) if i+1 == i1 else 1) > 68):
      s1 = s1[:-1] + '\n'
      col1 = 1

    s2 = s2 + ' '
    col2 = col2 + 1 + (len(str(n2-1)) if j == 0 else 1)
          
    if (col2 + (len(str(n2)) if j == 0 else 1) > 68):
      s2 = s2[:-1] + '\n'
      col2 = 1
            
    x = x+dx      
  
  y = y+dy
  
  # Prevent empty lines if the latest y-row ended at exactly the 68th column.
  if col1 != 1:
    s1 = s1 + '\n'
  if col2 != 1:
    s2 = s2 + '\n'
    
print s1
print "________________________________________________________________________"
print
print s2
