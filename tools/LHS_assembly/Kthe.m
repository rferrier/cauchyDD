function [K,C,ntot,node2c,c2node] = Kthe (nodes, elem, mat, order,...
    boundary, bc)
 % Computes the global thermic matrix
 % input : nodes    : list of nodes : [x1,y1;x2,y2;.....]
 %         elem     : list of elements
 %         mat      : Material : 0=isotropic, 1=orthotropic
 %         order    : order of the elements
 %         bc       : boundary conditions
 
 % mat = [0, E] if isotropic
 
 % output : K    : global stiffness matrix
 %          ntot : number of nodes on the uy=0 border

 % Only for debug purposes :
 %Kin = zeros(size(nodes,1), size(nodes,1));
 
 % First build the model's stiffness matrix
 mo = 0;
 mul = 1;
 
 if mat(1) == 0
    E   = mat(2);
    Sm1 = [E,0;0,E];
 else
    error('The required material behaviour is not implemented')
 end
 
 Kin = sparse(size(nodes,1), size(nodes,1)); % 'cause I love when things go fast
 
 for i=1:size(elem,1)
     Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
     nnodes = size(Xloc1,1);
     Xloc = zeros(nnodes,1);

     map = zeros(1,nnodes);
     for j=1:nnodes
         map(1,j) = elem(i,j);
         Xloc([2*j-1,2*j],1) = [Xloc1(j,1);Xloc1(j,2)];
     end
     
     Ke = thermat(Xloc,order,Sm1,0);
     
     % Build K
     Kin(map,map) = Kin(map,map) + Ke;
 end

 % First compute the number of elements with Dirichlet bc
 inodes = zeros(size(nodes,1),4); % Vector that stores witch coords are imposed
 ntot = 0;
 for i=1:size(boundary,1)
     entity = boundary(i,1);
     for j=1:size(bc,1)
         if bc(j,1) == entity
             for k=2:size(boundary,2)
                 inodes(boundary(i,k),1) = 1;
             end
         end
     end
 end
 
 for i=1:size(nodes,1)
     if inodes(i,1) ~= 0
         ntot = ntot+1;
     end
 end

 node2c = zeros( 1, size(nodes,1) ); % This list makes lines of C match with ddls
 c2node = zeros( 1, ntot );

 for j=1:size(bc,1)
     if bc(j,1) == 0 % Regularisation/x or y
         ntot = ntot+1;
     end
 end
 
 C = zeros(size(nodes,1),ntot);

 % Second time : populate C (I don't like growing matrixes)
 j=1;
 for i=1:size(nodes,1)
     if inodes(i,1) ~= 0
         C(i,j) = inodes(i,1);
         node2c(1,i) = j;
         c2node(1,j) = i;
         j=j+1;
     end
 end
 
 % Regularisation on x : sum(ux) = 0
 for k=1:size(bc,1)
     if bc(k,1) == 0 % Regularisation
         for i=1:size(nodes,1)
             C(i,j) = C(i,j) + 1;
         end
         j=j+1;
     end
 end

 K = [Kin,mul*C;mul*C',zeros(ntot)];
 
end
