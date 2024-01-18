function MCi =  MC_int(vtx_k,N)

%function to compute the integral of a function inside the k'th triangle
%with vertices at vtx_k using N sampoling points

MCi = zeros(N,size(vtx_k,2));

for i = 1:N
    r=rand; 
    s=rand;
    xa=1-sqrt(s); 
    xb=(1-r)*sqrt(s); 
    xc = 1-xa-xb;
    MCi(i,1) = xa*vtx_k(1,1) + xb*vtx_k(2,1) + xc*vtx_k(3,1);
    MCi(i,2) = xa*vtx_k(1,2) + xb*vtx_k(2,2) + xc*vtx_k(3,2);
end
