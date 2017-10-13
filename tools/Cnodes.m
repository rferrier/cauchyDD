function [ C, node2c, c2node, ntot ] = Cnodes ( nodes, index )
% This function computes the trace matrix C such that C'*u = u|index

 % First compute the number of elements with Dirichlet bc
 ntot = size(index,1);
 if ntot == 1
    ntot = size(index,2);
 end
 
 if size(index,1) > 1 && size(index,2) > 1
    warning('the index are not in a vector')
 end
 
 C = zeros(3*size(nodes,1),ntot);

 j=1;
 for i=1:ntot
    C(index(i),j) = 1;
    j = j+1;
 end

end
