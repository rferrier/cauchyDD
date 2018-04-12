function nodes2 = translateMesh( nodes1, x, y )
 % This function translates the given nocdes with a vector given by x and y
 
 nnodes = size(nodes1,1);
 Vector = [x*ones(nnodes,1), y*ones(nnodes,1)];
 nodes2 = nodes1 + Vector;

end

