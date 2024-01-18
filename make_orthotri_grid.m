function linmesh = make_orthotri_grid(n,L)

%Input: n the number of nodes in each direction over the inverval [0,1]
%
%Output: a grid definition consisted of orthogonal triangles in vertices
%and simplices matrices
%%
xp = linspace(0,L,n+1); yp = linspace(0,L,n+1);

[X,Y] = meshgrid(xp, yp);
vtx = [X(:),Y(:)];
DT = delaunayTriangulation(vtx);

vtx = DT.Points;
simp = DT.ConnectivityList;
%edg = edges(DT); %if needed
srf = freeBoundary(DT);
csimp = incenter(DT);

%Associate (interior) nodes to their corresponding elements
conn = vertexAttachments(DT);

%% Sort the rows of input
% simpl = sort(simpl,2);
% srfl = sort(srfl,2);
% 
% edges = [];
%  
% area = zeros(length(simpl),1);
% 
%  %loop over all linear elements
%  for i=1:size(simpl,1)
%      
%  %develop the 3 edges for each element
%  edges = [edges;[simpl(i,1), simpl(i,2),i]; [simpl(i,2),simpl(i,3),i];...
%          [simpl(i,1),simpl(i,3),i]];   
%      
%      
%  %compute the areas 
%  vtx_ith = vtxl(simpl(i,:),:); 
%  area(i) = polyarea(vtx_ith(:,1),vtx_ith(:,2));
%  end %for i
% % 
% % %Get the unique list of edges in the linear grid
%  [uedges,~,~] = unique(edges(:,1:2),'rows');
% % 
% % %Sort the rows of uedges
%  uedges = sort(uedges,2);

 
linmesh.vtx = vtx;
linmesh.simp = simp;
linmesh.srf = srf;
linmesh.csimp = csimp;
linmesh.area = triarea(vtx,simp);
linmesh.conn = conn;

