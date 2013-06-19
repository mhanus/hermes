import sys

if len(sys.argv) < 3 or len(sys.argv) > 5:
  print
  print "USAGE:"
  print
  print "  python gmsh2hermes.py INFILE OUTFILE [MAT_FILE] [BND_FILE]"
  print
  print "MAT_FILE : file with element marker names (one per line)"
  print "BND_FILE : file with boundary marker names (one per line)"
  print
  print "Integer physical groups as specified in GMSH will be used if"\
        " MAT_FILE or BND_FILE is not specified."
  print 
  sys.exit(1)
  
mat = dict()
bnd = dict()

vertices = []
elements = []
boundaries = []

glob_cnt = 0
if len(sys.argv) == 5:
  bndfile = open(sys.argv[4], 'r')

  for line in bndfile:
    s = line.split()
    if not s:
      continue
    
    if len(s) == 1:
      cnt = glob_cnt
      glob_cnt += 1
    else:
      cnt = int(s[0])
      
    bnd[cnt] = s[1]
    
  bndfile.close() 

glob_cnt = 0    
if len(sys.argv) >= 4:
  matfile = open(sys.argv[3], 'r')
  
  for line in matfile:
    s = line.split()
    if not s:
      continue
    
    if len(s) == 1:
      cnt = glob_cnt
      glob_cnt += 1
    else:
      cnt = int(s[0])
      
    mat[cnt] = s[1]
    
  matfile.close()

f = open(sys.argv[1], 'r')

for line in f:
  if '$MeshFormat' in line:
    state = 0
    next_state = 1
    continue
  elif '$Nodes' in line:
    state = 0
    next_state = 2
    continue
  elif '$Elements' in line:
    state = 0
    next_state = 4
    continue
    
  if state == 0:
    s = line.split()
    if not s:
      continue
    
    state = next_state
    
  if state == 1:
    s = line.split()
           
    if float(s[0]) != 2.2:
      print "Unsupported version of GMSH file. Anything may happen ;-)"
    
    state = 0
    next_state = 0 
    
  if state == 2:
    nvertices = int(line)
    state = 0
    next_state = 3
    
  if state == 3:
    if '$EndNodes' in line:
      continue
      
    s = line.split()
    
    conv = lambda x: float(x) if abs(float(x)) > 1e-11 else 0.0
    vertices.append(map(conv, s[1:-1]))
       
    state = 0
    next_state = 3
    
  if state == 4:
    nelems = int(line)
        
    state = 0
    next_state = 5 
  
  if state == 5:
    if '$EndElements' in line:
      continue    
    
    data = map(lambda x: int(x)-1, line.split())
        
    if data[1] == 0:    # 2-node line
      if len(bnd) > 0:
        ph_group = bnd[data[3]+1]
      else:
        ph_group = str(data[3]+1)
        
      boundaries.append(data[-2:]+[ph_group])
      
    elif data[1] == 1:  # 3-node triangle
      if len(mat) > 0:
        ph_group = mat[data[3]+1]
      else:
        ph_group = str(data[3]+1)
        
      elements.append(data[-3:]+[ph_group])
      
    elif data[1] == 2:  # 4-node quadrangle
      if len(mat) > 0:
        ph_group = mat[data[3]+1]
      else:
        ph_group = str(data[3]+1)
        
      elements.append(data[-4:]+[ph_group])
      
    else:
      print "Unsupported element type. Skipping ..."
    
    state = 0
    next_state = 5

f.close()

f = open(sys.argv[2], 'w')
stringize = lambda x: ("    " + str(x) + ",\n").replace("'", '"')

f.write("vertices = [\n")
f.writelines(map(stringize,vertices))
f.write("]\n\n")

f.write("elements = [\n")
f.writelines(map(stringize,elements))
f.write("]\n\n")

f.write("boundaries = [\n")
f.writelines(map(stringize,boundaries))
f.write("]\n")

f.close()
