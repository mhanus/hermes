clear all

pweights
num_elem = 7840;
N = length(pw)

k = 1;
I = zeros(num_elem*N,1);
J = zeros(num_elem*N,1);
V = zeros(num_elem*N,1);
for i=1:num_elem
  for j=1:N
    I(k) = i;
    J(k) = (j-1)*num_elem+i;
    V(k) = pw(j);
    k = k+1; 
  end
end

D = 8*pi*sparse(I, J, V);
save '-v7' 'D.mat' D

I = 1/(4*pi)*speye(num_elem);
M = repmat(I, N, 1);
save '-v7' 'M.mat' M
