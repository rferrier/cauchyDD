function f = lagr2forces2( lagr, c2node, entity, boundary, nnodes )
% This function reads the Lagrange multiplicators and deduces the reaction
% forces

% Inupt : lagr    : vector of Lagrange multiplicators
%         c2node  : Matrix of constraints ddl
%         entity  : The boundary on witch we need the forces (list)
%         nnodes  : The number of nodes

% Output : f     : Nodal reaction forces (list)

f = zeros(2*nnodes,size(entity,2));
for i=1:size(entity,2)
    enty = entity(1,i);
    fi = zeros(2*nnodes,1);
    fi(c2node,1) = -lagr(:,1) ;
    f(:,i) = keepField( fi, enty, boundary ); % Keep only on the desided boundary
end

end

