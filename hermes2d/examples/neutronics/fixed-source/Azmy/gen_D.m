clear all

pweights

ndof = 655360
N = length(pw)
ndof_per_ord = ndof/N

k = 1;
I = zeros(ndof,1);
J = zeros(ndof,1);
V = zeros(ndof,1);
for i=1:ndof_per_ord
  for j=1:N
    I(k) = i;
    J(k) = (j-1)*ndof_per_ord+i;
    V(k) = pw(j);
    k = k+1; 
  end
end

D = 8*pi*sparse(I, J, V);
save '-v7' 'D.mat' D

I = 1/(4*pi)*speye(ndof_per_ord);
M = repmat(I, N, 1);
save '-v7' 'M.mat' M
