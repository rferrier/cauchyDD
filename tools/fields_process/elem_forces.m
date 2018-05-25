function fe = elem_forces(f,nodes,boundary,entity,order)
% This function computes the forces per surface element from a Rhs on one boundary

 % Restrict boundary
 indices   = find(boundary(:,1)==entity);
 boundaryR = boundary(indices,:);

 nboun = size(boundaryR,1);
 nnodes = size(nodes,1);

 nf = size(f,2);

 if nboun == 0
    warning('Empty boundary');
 end

 if order ~= 1
   error('Only order 1 implemented');
 end

 if size(boundaryR,2) == 3  % 2D

    ntoelem = zeros(nnodes,1);
    for i=1:nboun
       ntoelem(boundaryR(i,2)) = ntoelem(boundaryR(i,2)) + 1;
       ntoelem(boundaryR(i,3)) = ntoelem(boundaryR(i,3)) + 1;
    end
    ntoelemm1 = zeros(nnodes,1);
    ntoelemm1(find(ntoelem)) = 1./ntoelem(find(ntoelem)); % Avoid the 0

    % Duplicate ntoelemm1
    ntoelemm = zeros(2*nnodes,1);
    ntoelemm(1:2:2*nnodes-1) = ntoelemm1;
    ntoelemm(2:2:2*nnodes)   = ntoelemm1;

    fR = f.*(ntoelemm*ones(1,nf));

    fe = zeros(2*size(boundary,1),nf);
    for i=1:nboun
       no1 = boundaryR(i,2); no2 = boundaryR(i,3);
       x1 = nodes(no1,1); y1 = nodes(no1,2);
       x2 = nodes(no2,1); y2 = nodes(no2,2);
       len = ((x2-x1)^2+(y2-y1)^2)^(.5);
       fe(2*indices(i)-1,:) = 1/len * (fR(2*no1-1,:) + fR(2*no2-1,:));
       fe(2*indices(i),:)   = 1/len * (fR(2*no1,:)   + fR(2*no2,:));
    end

 elseif size(boundaryR,2) == 4 % 3D
 
    ntoelem = zeros(nnodes,1);
    for i=1:nboun
       ntoelem(boundaryR(i,2)) = ntoelem(boundaryR(i,2)) + 1;
       ntoelem(boundaryR(i,3)) = ntoelem(boundaryR(i,3)) + 1;
       ntoelem(boundaryR(i,4)) = ntoelem(boundaryR(i,4)) + 1;
    end
    ntoelemm1 = zeros(nnodes,1);
    ntoelemm1(find(ntoelem)) = 1./ntoelem(find(ntoelem)); % Avoid the 0

    % Duplicate ntoelemm1
    ntoelemm = zeros(3*nnodes,1);
    ntoelemm(1:3:3*nnodes-2) = ntoelemm1;
    ntoelemm(1:3:3*nnodes-1) = ntoelemm1;
    ntoelemm(1:3:3*nnodes)   = ntoelemm1;

    fR = f.*(ntoelemm*ones(1,nf));

    fe = zeros(3*size(boundary,1),nf);
    for i=1:nboun
       no1 = boundaryR(i,2); no2 = boundaryR(i,3); no3 = boundaryR(i,4);
       x1 = nodes(no1,1); y1 = nodes(no1,2); z1 = nodes(no1,3);
       x2 = nodes(no2,1); y2 = nodes(no2,2); z2 = nodes(no2,3);
       x3 = nodes(no3,1); y3 = nodes(no3,2); z3 = nodes(no3,3);
       
       vpre = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1) ;...
                (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1) ;...
                (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) ];
       S = .5*abs(norm(vpre));
       
       fe(3*indices(i)-2,:) = 1/S * (fR(3*no1-2,:) + fR(3*no2-2,:) + fR(3*no3-2,:));
       fe(3*indices(i)-1,:) = 1/S * (fR(3*no1-1,:) + fR(3*no2-1,:) + fR(3*no3-1,:));
       fe(3*indices(i),:)   = 1/S * (fR(3*no1,:)   + fR(3*no2,:)   + fR(3*no3,:));
    end
    
 else
    error('Failed to identify the dimension');
 end

end
