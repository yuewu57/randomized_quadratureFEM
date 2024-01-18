function [areas] = triarea(vtx,simp)

for i=1:size(simp,1)
    vtx_i = [vtx(simp(i,:),:),ones(3,1)];
    areas(i) = abs(0.5*det(vtx_i));
end