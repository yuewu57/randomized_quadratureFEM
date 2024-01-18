function [outts] = dubs2(simp)

S=[];

for b=1:length(simp)
   x = simp(b,1);
   y = simp(b,2);
   z = simp(b,3);
   s1 = [x y];
   s2 = [x z];
   s3 = [y,z];
   Sn = [s1; s2; s3];
   S =  [S ; Sn];
end

N=sort(S,2);
M=sortrows(N);

count=1;
inc=1;

m=size(M,1);
for i=1:count:m-1
   
   ithrow = M(i,:);
   jthrow = M(i+1,:); 
   
   if ithrow == jthrow 
      M(i,:)=zeros(1,2);
      M(i+1,:)=zeros(1,2);
      count=2;
   end
   
   if ithrow ~= jthrow
   count=1;
   end
   
 end

M; 

for dk=1:m

if M(dk,:)~= zeros(1,2);
     XL(inc,:)= M(dk,:);
      inc=inc+1;
  end
   
end
XL;
outts = XL;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KURSK Ver 2 							  								%
% Developed by: Nick Polydorides 								%
% First year Ph.D. achievement									%
% Copyright (c) N. Polydorides, September 2000			   %	
% Required: MATLAB 5.3 or update									%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%