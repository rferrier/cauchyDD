function f = dirichletRhs( u, entity, C, boundary )
 % Builds the Dirichlet Rhght Hand Side making that u = ud on entity

 nnodes = size(C,1);
 f = zeros(size(C,1)+size(C,2), 1);
 
 for j=1:size(boundary,1)
    if boundary(j,1) == entity
        % Look for the column of C witch is not 0
        column = 0;
        for k=1:size(C,2)
            if C(2*boundary(j,2)-1,k) ~= 0
                column = k;
            end
        end
        if column ~= 0
            f( nnodes+column, 1 ) = u( 2*boundary(j,2)-1, 1 );
        end

        column = 0;
        for k=1:size(C,2)
            if C(2*boundary(j,2),k) ~= 0
                column = k;
            end
        end
        if column ~= 0
            f( nnodes+column, 1 ) = u( 2*boundary(j,2), 1 );
        end

        column = 0;
        for k=1:size(C,2)
            if C(2*boundary(j,3)-1,k) ~= 0
                column = k;
            end
        end
        if column ~= 0
            f( nnodes+column, 1 ) = u( 2*boundary(j,3)-1, 1 );
        end

        column = 0;
        for k=1:size(C,2)
            if C(2*boundary(j,3),k) ~= 0
                column = k;
            end
        end
        if column ~= 0
            f( nnodes+column, 1 ) = u( 2*boundary(j,3), 1 );
        end
    end
 end
end

