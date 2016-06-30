function f = dirichletRhs2( u, entity, c2node, boundary, nnodes, varargin )
 % Builds the Dirichlet Rhght Hand Side making that u = ud on entity
 
 % component to keep
 co = 0;
 if numel(varargin)>0
     co = cell2mat(varargin(1));
 end

 up = keepField( u, entity, boundary, co ); % Keep only on the desided boundary
 f = [ zeros(2*nnodes,1) ; up(c2node, 1) ];
 
end

