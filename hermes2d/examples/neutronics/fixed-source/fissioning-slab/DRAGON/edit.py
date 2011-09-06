# Number of mesh subdivisions in each of the two regions.
nx = 50
ny = 50
ax = 5.
ay = 5.

# Length of subdivided mesh intervals.
dx = ax / nx
dy = ay / ny

# Resulting strings with edit regions.
s = ''
# Gradually increased number representing merge regions in the tally.
n = 1
# Column counter for wrapping after 68th column.
col = 1

y = 0.
for j in range(0,ny):

  x = 0.
  for i in range(0,nx):
    
    # Diagonal flux.  
    if (i == j):
      s = s + str(n)
      n = n + 1
      print str(x+dx/2.) + ", " + str(y+dy/2.)
    else:
      s = s + '0'
      
    # Append a space and either the merge-region index or zero (indicating that
    # the current mesh region is to be ignored in this tally).
    s = s + ' '
    col = col + 1 + (len(str(n-1)) if i == j else 1)
    
    # Should the next entry overflow the 68th column, remove the trailing space
    # and break the line.
    if (col + (len(str(n)) if i+1 == j else 1) > 68):
      s = s[:-1] + '\n'
      col = 1
            
    x = x+dx      
  
  y = y+dy
  
  # Prevent empty lines if the latest y-row ended at exactly the 68th column.
  if col != 1:
    s = s + '\n'
    
print s
