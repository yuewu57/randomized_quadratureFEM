%Set the size of the square 
L=1;

%Set the number of elements on L
n=2^6;

h=L/n;


%Generate mesh
mesh = make_orthotri_grid(n,L);

%Make gradients matrix D and volumes diagonal Z
[D,Z] = make_laplacian(mesh.vtx,mesh.simp);

%Extract boundary nodes to impose the boundary condition
nbd = unique(mesh.srf(:));
ndof = setdiff(1:length(mesh.vtx),nbd);

%Apply Dirichlet boundary conditions
Dc = D(:,ndof);

%The Laplacian matrix
%A = Dc.'*Z*Dc;

clear D 


%Set the sigma function sigma(x,y) = abs(sin(x*y))+0.2;
sigma_MC = zeros(size(mesh.simp,1),1);

for k=1:size(mesh.simp,1)
    %100 integration points uniformly distributed on the k'th simplex
    MCi =  MC_int(mesh.vtx(mesh.simp(k,:),:),1);
    sigma_MC(k) = sigma(MCi);
    
end

Z_MC = Z.*sparse(1:2*size(mesh.simp,1),1:2*size(mesh.simp,1),kron(sigma_MC,ones(2,1)));

A_MC = Dc.'*Z_MC*Dc;


sigma_MC_c = zeros(size(mesh.simp,1),1);

for k=1:size(mesh.simp,1)
    %100 integration points uniformly distributed on the k'th simplex
    MCi_c =  MC_int(mesh.vtx(mesh.simp(k,:),:),100);
    sigma_MC_c(k) = sigma(MCi_c);
    
end

Z_MC_c = Z.*sparse(1:2*size(mesh.simp,1),1:2*size(mesh.simp,1),kron(sigma_MC_c,ones(2,1)));

A_MC_c = Dc.'*Z_MC_c*Dc;
% trace the element connectivity of interior nodes (node-2-element) map
n2e = cell2mat(mesh.conn(ndof));


M=50;

%the matrix with the MC simulated rhs vectors
F = zeros(length(ndof),M);

for m=1:M %repeat M times


for i=1:length(ndof)

%node index for the i'th dof
i_dof = ndof(i);

%element indices including the i_dof
i_ele = n2e(i,:);

%initialise the 6 points matrix
Z_i = zeros(6,2);
Fun_Z_i=zeros(6,1);

for j=1:length(i_ele) %this loop is for computing (Z_i_ele)_j f

%the area of the jth element
i_ele_area = mesh.area(i_ele(j));

%the nodes in the i_ele
nodes_i_ele = mesh.simp(i_ele(j),:);

%find the index of the i_dof in the element definition
i_dof_ind = find(nodes_i_ele==i_dof);

%The coordinates of the i_ele(j) element
j_ele_coords = mesh.vtx(nodes_i_ele,:);

%The map from global to local coordinates. Coeffs for constant, alpha
%and beta
x_a_b = [j_ele_coords(1,1),j_ele_coords(2,1)-j_ele_coords(1,1),j_ele_coords(3,1)-j_ele_coords(1,1)];
y_a_b = [j_ele_coords(1,2),j_ele_coords(2,2)-j_ele_coords(1,2),j_ele_coords(3,2)-j_ele_coords(1,2)];
%end %for j elements connected to i_dof

%Select the appropriate local function definition on the standard
%simplex. Develop the 3 integrands
switch i_dof_ind
case 1 %\hat \phi(alpha,beta) = 1 - alpha - beta
%% Sampling process
z1 = rand(1); z2 = rand(1);
%Check if (z1,z2) are within the standard simplex
if z1+z2 <= 1
Z_i_dof_j = [z1,z2];
else
Z_i_dof_j = [1-z1,1-z2];
end
case 2 %\hat \phi(alpha,beta) = alpha
%% Sampling process
z1 = rand(1); z2 = rand(1);
%Check if (z1,z2) are within the standard simplex
if z1+z2 <= 1
Z_i_dof_j = [z1,z2];
else
Z_i_dof_j = [1-z1,1-z2];
end
case 3 %\hat \phi(alpha,beta) = beta
%% Sampling process
z1 = rand(1); z2 = rand(1);
%Check if (z1,z2) are within the standard simplex
if z1+z2 <= 1
Z_i_dof_j = [z1,z2];
else
Z_i_dof_j = [1-z1,1-z2];
end
end
%the Z_i^{j} sample in global coordinates is
Z_i(j,:) = [dot(x_a_b,[1,Z_i_dof_j(1),Z_i_dof_j(2)]),dot(y_a_b,[1,Z_i_dof_j(1),Z_i_dof_j(2)])];
Fun_Z_i(j)=i_ele_area*Fun1(Z_i(j,1) ,Z_i(j,2));
end %j=6


F(ndof(i),m) = sum(Fun_Z_i);
%figure; trimesh(mesh.simp(i_ele,:),mesh.vtx(:,1),mesh.vtx(:,2));
%hold on; plot(Z_i(:,1),Z_i(:,2),'x');
%plot(Z_if(:,1),Z_if(:,2),'ro');
end %i all ndofs

end %m=M

%Extract F at the degrees of freedom - interior nodes
Fn = F(ndof,:);
%Solve the model fem eq
U_MC=A_MC\Fn;

U_MC_c=A_MC_c\Fn;

%Test

u1 = zeros(length(mesh.vtx),1);
u1(ndof) = U_MC(:,3);
figure; trisurf(mesh.simp,mesh.vtx(:,1),mesh.vtx(:,2),u1);

u2=zeros(length(mesh.vtx),1);
u2(ndof) = mean(U_MC,2);
figure; trisurf(mesh.simp,mesh.vtx(:,1),mesh.vtx(:,2),u2);

u3=zeros(length(mesh.vtx),1);
u3(ndof) = mean(U_MC_c,2);
figure; trisurf(mesh.simp,mesh.vtx(:,1),mesh.vtx(:,2),u3);
