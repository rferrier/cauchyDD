function f = lagr2forces( lagr, C, entity, boundary )
% This function reads the Lagrange multiplicators and deduces the reaction
% forces

% Inupt : lagr   : vector of Lagrange multiplicators
%         C      : Matrix of constraints ddl
%         entity : The boundary on witch we need the forces (list)

% Output : f     : Nodal reaction forces (list)

f = zeros(size(C,1),size(entity,2));
for i=1:size(entity,2)
    enty = entity(1,i);
    for j=1:size(boundary,1)
        if boundary(j,1) == enty
            % Look for the column of C witch is not 0
            column = 0;
            for k=1:size(C,2)
                if C(2*boundary(j,2)-1,k) ~= 0
                    column = k;
                end
            end
            if column ~= 0
                f( 2*boundary(j,2)-1, i ) = - lagr(column,1);
            end
            
            column = 0;
            for k=1:size(C,2)
                if C(2*boundary(j,2),k) ~= 0
                    column = k;
                end
            end
            if column ~= 0
                f( 2*boundary(j,2), i ) = - lagr(column,1);
            end
            
            column = 0;
            for k=1:size(C,2)
                if C(2*boundary(j,3)-1,k) ~= 0
                    column = k;
                end
            end
            if column ~= 0
                f( 2*boundary(j,3)-1, i ) = - lagr(column,1);
            end
            
            column = 0;
            for k=1:size(C,2)
                if C(2*boundary(j,3),k) ~= 0
                    column = k;
                end
            end
            if column ~= 0
                f( 2*boundary(j,3), i ) = - lagr(column,1);
            end
        end
    end
end

end

