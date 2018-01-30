function M = rigidModes(nodes)
% This function computes the 2D mecanical rigid modes associated with a mesh

nnodes = size(nodes,1); i = 1:nnodes;

Mx = zeros(2*nnodes,1); My = zeros(2*nnodes,1);
Mx(1:2:2*nnodes-1) = 1; My(2:2:2*nnodes) = 1;

Mtz = zeros(2*nnodes,1);
Mtz(2*i-1) = -nodes(:,2); Mtz(2*i) = nodes(:,1);

M = [Mx,My,Mtz];
end

