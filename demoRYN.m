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
A = Dc.'*Z*Dc;

clear D Dc Z
% trace the element connectivity of interior nodes (node-2-element) map
n2e = cell2mat(mesh.conn(ndof));

%We assume a choice of f(x,y) on [0,L]x[0,L] as


M=30;

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
            SS=1; %enter sampling process
                while SS > 0
   
                    z1 = rand(1); z2 = rand(1);
                    %Check if (z1,z2) are within the standard simplex
                    if z1+z2 <= 1
                        Z_i_dof_j = [z1,z2];
                    else    
                        Z_i_dof_j = [1-z1,1-z2];
                    end
                    C = 3*rand(1);    %small c = 3
                    if (2*C - 6*(1 - Z_i_dof_j(1) - Z_i_dof_j(2)))<= 0 %sample valid
                    SS=0;
                    end
                end %while sampling process
            
    
        case 2 %\hat \phi(alpha,beta) = alpha
             %% Sampling process
             SS=1; %enter sampling process
                while SS > 0
                    z1 = rand(1); z2 = rand(1);
                    %Check if (z1,z2) are within the standard simplex
                    if z1+z2 <= 1
                       Z_i_dof_j = [z1,z2];
                    else    
                       Z_i_dof_j = [1-z1,1-z2];
                    end
                    C = 3*rand(1);    %small c = 3
                    if (2*C - 6*(Z_i_dof_j(1)))<= 0 %sample valid
                    SS=0;
                    end
                end %while sampling process 
        case 3 %\hat \phi(alpha,beta) = beta
        
            %% Sampling process
            SS=1; %enter sampling process
                while SS > 0
                    z1 = rand(1); z2 = rand(1);
                    %Check if (z1,z2) are within the standard simplex
                    if z1+z2 <= 1
                        Z_i_dof_j = [z1,z2];
                    else    
                        Z_i_dof_j = [1-z1,1-z2];
                    end
                    C = 3*rand(1);    %small c = 3
                    if (2*C - 6*(Z_i_dof_j(2)))<= 0 %sample valid
                    SS=0;
                    end
                end %while sampling process 
        
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
U = A\Fn;

%Test
u = zeros(length(mesh.vtx),1);
u(ndof) = U(:,3);
figure; trisurf(mesh.simp,mesh.vtx(:,1),mesh.vtx(:,2),u);


u2=zeros(length(mesh.vtx),1);
u2(ndof) = mean(U,2);
figure; trisurf(mesh.simp,mesh.vtx(:,1),mesh.vtx(:,2),u2);

q=0.49;

%One point Gaussian quadrature rule for %f(x,y) = sin(2*pi*x)+
%abs(cos(2*pi*y))+
%f(x,y) =|x-y|^(-q)+10ysign(2y-x). This needs to be multiplied with \phi(x)
f1q = zeros(size(mesh.vtx,1),1);
for k=1:size(mesh.simp,1)
        V = mesh.vtx(mesh.simp(k,:),:);          
        f_a_b_c = abs(((V(3,2)-V(1,2))/3 + (V(2,2)-V(1,2))/3 + V(1,2))...
            -((V(3,1)-V(1,1))/3 + (V(2,1)-V(1,1))/3 + V(1,1))+eps).^(-q)...
            +10*((V(3,2)-V(1,2))/3 + (V(2,2)-V(1,2))/3 + V(1,2))...
          *sign(2*((V(3,2)-V(1,2))/3 + (V(2,2)-V(1,2))/3 + V(1,2))...
          -((V(3,1)-V(1,1))/3 + (V(2,1)-V(1,1))/3 + V(1,1))); 
        %sin(2*pi*((V(3,1)-V(1,1))/3 + (V(2,1)-V(1,1))/3 + V(1,1))) ...
         % + abs(cos(2*pi*((V(3,2)-V(1,2))/3 + (V(2,2)-V(1,2))/3 + V(1,2))))...
   
        %1-point guass quadrature rule
        f1q(mesh.simp(k,:)) = f1q(mesh.simp(k,:)) + (mesh.area(k)/3)*f_a_b_c;
end
 
%Check the f1q to see if it makes sense
figure; trisurf(mesh.simp,mesh.vtx(:,1),mesh.vtx(:,2),f1q);

%Extract the rhs at the interior nodes
f1qndof = f1q(ndof); mF = mean(Fn,2);
figure; plot(mF,'.r'); hold on; plot(f1qndof,'.');

U1=A\f1qndof;
u1 = zeros(length(mesh.vtx),1);
u1(ndof) = U1;
figure; trisurf(mesh.simp,mesh.vtx(:,1),mesh.vtx(:,2),u1);

%Compute the mass matrix - bilinear form (phi_i,phi_j) - each local matrix
%is (1/6)*area at the diagonal and (1/12)*area off the diagonal
B = spalloc(length(mesh.vtx),length(mesh.vtx), 9*length(mesh.vtx));
for i=1:size(mesh.simp,1)
    B(mesh.simp(i,:),mesh.simp(i,:)) = B(mesh.simp(i,:),mesh.simp(i,:))...
        +mesh.area(i)*[1/6, 1/12, 1/12; 1/12, 1/6, 1/12; 1/12, 1/12, 1/6];
end

%Extract the mass matrix at the dofs 
Bx = B(ndof,ndof);
