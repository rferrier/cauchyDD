function [K,C,ntot,node2c,c2node] = KrigDam (nodes, elem, E, nu, order,...
    boundary, bc, u, d)
 % Computes the global tangent stiffness matrix for damage modes
 % input : nodes    : list of nodes : [x1,y1;x2,y2;.....]
 %         elem     : list of elements
 %         E        : Young Modulus
 %         nu       : Poisson ratio
 %         order    : order of the elements
 %         bc       : boundary conditions
 %         u        : previous displacment
 %         d        : damage value
 
 % output : K    : global stiffness matrix
 %          ntot : number of nodes on the uy=0 border

 % Only for debug purposes :
 %Kin = zeros(2*size(nodes,1), 2*size(nodes,1));
 
 % First build the model's stiffness matrix
 S = 1/E*[1,-nu,0 ; -nu,1,0 ; 0,0,2*(1+nu)];
 Sm1 = inv(S);
 
 Kin = sparse(2*size(nodes,1), 2*size(nodes,1)); % 'cause I love when things go fast
 
 for i=1:size(elem,1)
     Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
     nnodes = size(Xloc1,1);
     Xloc = zeros(2*nnodes,1);
     uloc = zeros(2*nnodes,1);      % Extract local displacement
     for j=1:nnodes
         Xloc([2*j-1,2*j],1) = [Xloc1(j,1);Xloc1(j,2)];
         uloc([2*j-1,2*j],1) = [ 2*u(elem(i))-1 ; u(elem(i)) ];
     end
     
     % Compute tr(epsilon) : TODO : compute epsilon 33 
     % ('cause right now, it's false)
     eps1 = 1/3*( (uloc(3)-uloc(1))/(Xloc(3)-Xloc(1)) +...
         (uloc(5)-uloc(1))/(Xloc(5)-Xloc(1)) +...
         (uloc(5)-uloc(3))/(Xloc(5)-Xloc(3)) );
     eps2 = 1/3*( (uloc(4)-uloc(2))/(Xloc(4)-Xloc(2)) +...
         (uloc(6)-uloc(2))/(Xloc(6)-Xloc(2)) +...
         (uloc(6)-uloc(4))/(Xloc(6)-Xloc(4)) );
     
     treps = 1/3*(eps1+eps2);
     
     if treps >= 0 % traction
        Ke = stifmat(Xloc,order,Sm1,0);
     else % compression
        Ke = stifmat(Xloc,order,Sm1.*(1-d),0);
     end
     
     % Build K
     map = [2*elem(i,1)-1,2*elem(i,1),2*elem(i,2)-1,...
            2*elem(i,2),2*elem(i,3)-1,2*elem(i,3)];
     Kin(map,map) = Kin(map,map) + Ke;
 end
 
%  % Boundary conditions with Lagrange multiplicators
%  d_lines = [];
%  for i = 1:size(bc,1)
%      d_lines(i,1) = bc(i,1);
%  end

 % First compute the number of elements with Dirichlet bc
 inodes = zeros(size(nodes,1),4); % Vector that stores witch coords are imposed
 ntot = 0;
 for i=1:size(boundary,1)
     entity = boundary(i,1);
     for j=1:size(bc,1)
         if bc(j,1) == entity
             for k=2:size(boundary,2)

                 if bc(j,2) == 1     % Mark this node as x
                     
%                     if inodes(boundary(i,k),3) ~= 0 &&...
%                             inodes(boundary(i,k),3) ~= entity
%                         inodes(boundary(i,k),1) = inodes(boundary(i,k),1) + 1;
%                     else
                        inodes(boundary(i,k),1) = 1;
%                     end
%                     inodes(boundary(i,k),3) = entity;
                 elseif bc(j,2) == 2 % Mark this node as y
                     
%                     if inodes(boundary(i,k),4) ~= 0 &&...
%                             inodes(boundary(i,k),4) ~= entity
%                         inodes(boundary(i,k),2) = inodes(boundary(i,k),2) + 1;
%                     else
                        inodes(boundary(i,k),2) = 1;
%                     end
%                     inodes(boundary(i,k),4) = entity;
                 else
                     error('Unable to read Dirichlet BC axis')
                 end
                 
             end
         end
     end
 end
 
%  % Detect if a node belongs to 2 Dirichlet's BC
%  for i=1:size(bc,1)
%  end
 
 for i=1:size(nodes,1)
     if inodes(i,1) ~= 0
         ntot = ntot+1;
     end
     if inodes(i,2) ~= 0
         ntot = ntot+1;
     end
 end

 node2c = zeros( 1, 2*nnodes ); % This list makes lines of C match with ddls
 c2node = zeros( 1, ntot );

 for j=1:size(bc,1)
     if bc(j,1) == 0 % Regularisation/x or y
         ntot = ntot+1;
     end
 end
 
 C = zeros(2*size(nodes,1),ntot);

 % Second time : populate C (I don't like growing matrixes)
 j=1;
 for i=1:size(nodes,1)
     if inodes(i,1) ~= 0
         C(2*i-1,j) = inodes(i,1);
         node2c(1,2*i-1) = j;
         c2node(1,j) = 2*i-1;
         j=j+1;
     end
     if inodes(i,2) ~= 0
         C(2*i,j) = inodes(i,2);
         node2c(1,2*i) = j;
         c2node(1,j) = 2*i;
         j=j+1;
     end
 end
 
 % Regularisation on x : sum(ux) = 0
 for k=1:size(bc,1)
     if bc(k,1) == 0 && bc(k,2) == 1 % Regularisation/x
         for i=1:size(nodes,1)
             C(2*i-1,j) = C(2*i-1,j) + 1;
         end
         j=j+1;
     end
     if bc(k,1) == 0 && bc(k,2) == 2 % Regularisation/y
         for i=1:size(nodes,1)
             C(2*i,j) = C(2*i,j) + 1;
         end
         j=j+1;
     end
     if bc(k,1) == 0 && bc(k,2) == 3 % Regularisation/theta
         % First, find the ~=barycenter of the solid
         x0 = sum( nodes(:,1) ) / size(nodes,1);
         y0 = sum( nodes(:,2) ) / size(nodes,1);
         for i=1:size(nodes,1)
             x = nodes(i,1);
             y = nodes(i,2);
             if abs(x-x0) > 1e-1 % avoid 1/small
                 C(2*i,j) = C(2*i,j) + 1/(x-x0);
             end
             if abs(y-y0) > 1e-1
                 C(2*i-1,j) = C(2*i-1,j) - 1/(y-y0);
             end
             
         end
         j=j+1;
     end
 end

 K = [Kin,C;C',zeros(ntot)];
 
end