function [ nodes,elements,ntoelem,boundary,order ] =...
    createmesh( X, Y, dx, dy, varargin )
% This function cretates the regular mesh corresponding to a rectangle
 nx = floor(X/dx)+1;
 ny = floor(Y/dy)+1;
 dx = X/(nx-1);
 dy = Y/(ny-1);

 nnodes = nx*ny + (nx-1)*(ny-1);
 nelem  = 4*(nx-1)*(ny-1);
 nbound = 2*(nx+ny-2);
 
 nodes = zeros(nnodes,2);
 elements = zeros(nelem,3);
 boundary = zeros(nbound,3);
 order    = 1;

 % First : the nodes at the vertex
 for i=1:ny
     for j=1:nx
         nodes(j+(i-1)*nx,:) = [(j-1)*dx, (i-1)*dy];
     end
 end

 % Second : the nodes in center
 for i=1:ny-1
     for j=1:nx-1
         nodes(nx*ny+j+(i-1)*(nx-1),:) = [(j-1)*dx+dx/2, (i-1)*dy+dy/2];
     end
 end

 % The elements
 for i=1:ny-1
     for j=1:nx-1
         elements(4*(j-1)+4*(i-1)*(nx-1)+1,:) =...
             [j+(i-1)*nx, j+1+(i-1)*nx, nx*ny+j+(i-1)*(nx-1)];
         elements(4*(j-1)+4*(i-1)*(nx-1)+2,:) =...
             [j+(i-1)*nx, j+(i)*nx, nx*ny+j+(i-1)*(nx-1)];
         elements(4*(j-1)+4*(i-1)*(nx-1)+3,:) =...
             [j+1+(i)*nx, j+1+(i-1)*nx, nx*ny+j+(i-1)*(nx-1)];
         elements(4*(j-1)+4*(i-1)*(nx-1)+4,:) =...
             [j+(i)*nx, j+1+(i)*nx, nx*ny+j+(i-1)*(nx-1)];
     end
 end

 % ntoelem is a table storing the number of elements to witch each node belongs
 ntoelem = zeros(size(nodes,1), 1);
 for i=1:size(elements,1)
     for j=1:3
         ntoelem(elements(i,j),1) = ntoelem(elements(i,j),1) + 1;
     end
 end
 
 % The boundary elements
 for i=1:nx-1
     boundary(i,:) = [1,i,i+1];
     boundary(i+ny+nx-2,:) = [3,(ny-1)*nx+i,(ny-1)*nx+i+1];
 end
 for i=1:ny-1
     boundary(i+nx-1,:) = [2,i*nx,(i+1)*nx];
     boundary(i+ny+2*nx-3,:) = [4,(i-1)*(nx)+1,(i)*(nx)+1];
 end
 
 % Display the mesh in a .msh file (for GMSH)
 if numel( varargin ) > 0
     name = cell2mat( varargin(1) );
     
    fmid = fopen(name,'w');
    fprintf(fmid,'%s\n%s\n%s\n','$MeshFormat','2.2 0 8','$EndMeshFormat');
     
    % Nodes
    fprintf(fmid,'%s\n','$Nodes');
    fprintf(fmid,'%d\n',nnodes);
    for n=1:nnodes
        fprintf(fmid,'%d %d %d %d \n',n,nodes(n,1),nodes(n,2),0);   
    end
    fprintf(fmid,'%s\n','$EndNodes');
    
    % Elements
    fprintf(fmid,'%s\n','$Elements');
    fprintf(fmid,'%d\n',nelem+nbound);
    for n=1:nbound
        fprintf(fmid,'%d %d %d %d %d %d %d \n',n,1,2,boundary(n,1),...
            boundary(n,1),boundary(n,2),boundary(n,3));   
    end
    for n=1:nelem
        fprintf(fmid,'%d %d %d %d %d %d %d %d \n',n+nbound,2,2,5,1,...
            elements(n,1),elements(n,2),elements(n,3));   
    end
    fprintf(fmid,'%s\n','$EndElements');
 end
 
end
