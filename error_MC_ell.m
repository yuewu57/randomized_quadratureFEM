%Set the size of the square 
L=1;

%Set the number of elements on L
n=2^8;

h=L/n;

q=0.49;


%Generate mesh
mesh = make_orthotri_grid(n,L);

%Make gradients matrix D and volumes diagonal Z
[D,Z] = make_laplacian(mesh.vtx,mesh.simp);

%Extract boundary nodes to impose the boundary condition
nbd = unique(mesh.srf(:));
ndof = setdiff(1:length(mesh.vtx),nbd);

%Apply Dirichlet boundary conditions
Dc = D(:,ndof);


clear D 


%Set the sigma function;
sigma_MC = zeros(size(mesh.simp,1),1);

for k=1:size(mesh.simp,1)
    MCi =  MC_int(mesh.vtx(mesh.simp(k,:),:),1);
    sigma_MC(k) = sigma(MCi);
    
end

Z_MC = Z.*sparse(1:2*size(mesh.simp,1),1:2*size(mesh.simp,1),kron(sigma_MC,ones(2,1)));

A_MC = Dc.'*Z_MC*Dc;


sigma_MC_bench = zeros(size(mesh.simp,1),1);
for k=1:size(mesh.simp,1)
    %100 integration points uniformly distributed on the k'th simplex
    MCi_bench =  MC_int(mesh.vtx(mesh.simp(k,:),:),1000);
    sigma_MC_bench(k) = sigma(MCi_bench);
    
end

Z_MC_bench = Z.*sparse(1:2*size(mesh.simp,1),1:2*size(mesh.simp,1),kron(sigma_MC_bench,ones(2,1)));

A_MC_bench = Dc.'*Z_MC_bench*Dc;
% trace the element connectivity of interior nodes (node-2-element) map
n2e = cell2mat(mesh.conn(ndof));

%We assume a choice of f(x,y) on [0,L]x[0,L] as
%f(x,y) = sin(cs*x)+ abs(cos(cc*y)) for some constants cs and cc
%cs = 2^9*pi; cc = 2^5*pi;

M=3000;

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
    Fun_Z_i(j)=i_ele_area*Fun_rough(Z_i(j,1) ,Z_i(j,2));
    end %j=6
    

%figure; trimesh(mesh.simp(i_ele,:),mesh.vtx(:,1),mesh.vtx(:,2));
%hold on; plot(Z_i(:,1),Z_i(:,2),'x');
%plot(Z_if(:,1),Z_if(:,2),'ro');
F(ndof(i),m) = sum(Fun_Z_i);


end %i all ndofs

end %m=M

Fn = F(ndof,:);
%Solve the model fem eq
U = A_MC\Fn;

%error term in strong sense, i.e., expectation of error in H_0^1(D) norm
error1 = 0;
error2 = 0;

for k=1:M
 error1 = error1+ U(:,k)'*A_MC_bench*U(:,k);
 if k<M
 for kk=k+1:M
     error2=error2 +U(:,k)'*A_MC_bench*U(:,kk);
 end
 end
end
norm_error=(error1-error2*2/M)/(M-1);

sqrt(norm_error)

%Compute the mass matrix - bilinear form (phi_i,phi_j) - each local matrix
%is (1/6)*area at the diagonal and (1/12)*area off the diagonal
%B = spalloc(length(mesh.vtx),length(mesh.vtx), 9*length(mesh.vtx));
%for i=1:size(mesh.simp,1)
 %   B(mesh.simp(i,:),mesh.simp(i,:)) = B(mesh.simp(i,:),mesh.simp(i,:))...
  %      +mesh.area(i)*[1/6, 1/12, 1/12; 1/12, 1/6, 1/12; 1/12, 1/12, 1/6];
%end

%Extract the mass matrix at the dofs 
%Bx = B(ndof,ndof);

%Compute the error in L2 norm 
%error3 = 0;
%error4 = 0;

%for k=1:M
 %error3 = error3+ U(:,k)'*Bx*U(:,k);
%  if k<M
%  for kk=k+1:M
%      error4=error4 +U(:,k)'*Bx*U(:,kk);
%  end
%  end
% end
%L_error=(error3-error4*2/M)/(M-1);

%sqrt(L_error)
