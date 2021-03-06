function f = dirichletRhs2( u, entity, c2node, boundary, nnodes, varargin )
 % Builds the Dirichlet Rhght Hand Side making that u = ud on entity
 
 % component to keep
 co = 0;
 if numel(varargin)>0
     co = cell2mat(varargin(1));
 end

 mul = 1;
 if numel(varargin)>1
     mul = cell2mat(varargin(2));
 end
 
 s2 = size(u,2);
 
 up = keepField( u, entity, boundary, co ); % Keep only on the desided boundary
 f = [ zeros(2*nnodes,s2) ; mul*up(c2node, :) ];
 
end

