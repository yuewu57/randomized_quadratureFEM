function [D,Z] = make_laplacian(vtx,simp)

%Based on gm_assemble by S.Vavasis , http://www.cs.cornel.edu/ 

[vr, vc] = size(vtx);
[sr, sc] = size(simp);
a = ones(length(simp),1);

ilist = kron((1:vc*sr), [1,1]);
jlist = zeros(1,sr*vc*2);
slist = zeros(1,sr*vc*2);

for d = 1 : vc
  jlist(2 * d - 1 : 2 * vc : sr * vc * 2) = simp(:,1)';
  jlist(2 * d : 2 * vc : sr * vc * 2) = simp(:, d + 1);
  slist(2 * d - 1 : 2 * vc : sr * vc * 2) = -ones(1,sr);
  slist(2 * d : 2 * vc : sr * vc * 2) = ones(1,sr);
end

A0 = sparse(ilist,jlist,slist,vc*sr,vr);


if vc == 2
  J1 = A0 * vtx(:,1);
  J2 = A0 * vtx(:,2);
  J11 = J1(1:2:sr*2);
  J12 = J1(2:2:sr*2);
  J21 = J2(1:2:sr*2);
  J22 = J2(2:2:sr*2);
  detJ = J11 .* J22 - J21 .* J12;
  
  
  invJ11 = J22 ./ detJ;
  invJ12 = -J12 ./ detJ;
  invJ21 = -J21 ./ detJ;
  invJ22 = J11 ./ detJ;
elseif vc == 3
  J1 = A0 * vtx(:,1);
  J2 = A0 * vtx(:,2);
  J3 = A0 * vtx(:,3);
  J11 = J1(1:3:sr*3);
  J12 = J1(2:3:sr*3);
  J13 = J1(3:3:sr*3);
  J21 = J2(1:3:sr*3);
  J22 = J2(2:3:sr*3);
  J23 = J2(3:3:sr*3);
  J31 = J3(1:3:sr*3);
  J32 = J3(2:3:sr*3);
  J33 = J3(3:3:sr*3);
  detJ = J11 .* J22 .* J33 - J11 .* J23 .* J32 - J12 .* J21 .* J33 ...
          + J12 .* J23 .* J31 + J13 .* J21 .* J32 - J13 .* J22 .* J31;
       
       
  invJ11 = (J22 .* J33 - J23 .* J32) ./ detJ;
  invJ12 = (J32 .* J13 - J12 .* J33) ./ detJ;
  invJ13 = (J12 .* J23 - J22 .* J13) ./ detJ;
  invJ21 = (J31 .* J23 - J21 .* J33) ./ detJ;
  invJ22 = (J11 .* J33 - J13 .* J31) ./ detJ;
  invJ23 = (J21 .* J13 - J11 .* J23) ./ detJ;
  invJ31 = (J21 .* J32 - J31 .* J22) ./ detJ;
  invJ32 = (J31 .* J12 - J11 .* J32) ./ detJ;
  invJ33 = (J11 .* J22 - J21 .* J12) ./ detJ;
else
  error('Master matrix construction failed')
end


ilist = kron((1 : vc * sr), ones(1,vc));
jlist = zeros(1, sr*vc^2);
for d = 1 : vc 
  jlist(d:vc:sr*vc^2) = kron((d:vc:vc*sr),ones(1,vc));
end

if vc == 2
  slist = zeros(1,sr*4);
  slist(1:4:sr*4) = invJ11;
  slist(2:4:sr*4) = invJ21;
  slist(3:4:sr*4) = invJ12;
  slist(4:4:sr*4) = invJ22;
else
  slist = zeros(1,sr*9);
  slist(1:9:sr*9) = invJ11;
  slist(2:9:sr*9) = invJ21;
  slist(3:9:sr*9) = invJ31;
  slist(4:9:sr*9) = invJ12;
  slist(5:9:sr*9) = invJ22;
  slist(6:9:sr*9) = invJ32;
  slist(7:9:sr*9) = invJ13;
  slist(8:9:sr*9) = invJ23;
  slist(9:9:sr*9) = invJ33;
end


ElJac = sparse(ilist,jlist,slist,vc*sr,vc*sr);
D = ElJac * A0;


Vols = abs(detJ) / prod(1:vc);

materials = length(a);
volumes = size(Vols);

if materials ~= volumes
  error('Some elements have not been assigned');
end

Z = sparse( (1:vc*sr), (1:vc*sr),kron( (a .* Vols).',ones(1,vc)) );

%Ef = D'*Z*D; 

 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KURSK Ver 2 							  								%
% Developed by: Nick Polydorides 								%
% First year Ph.D. achievement									%
% N. Polydorides, September 2000			   %	
% Required: MATLAB 5.3 or update									%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%