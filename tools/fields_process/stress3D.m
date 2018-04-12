function sigma = stress3D( u,mat,nodes,elem,order,gorp,ntoelem )
 % This function computes the stresses in the material
 % input : u        : global solution in displacement
 %         E        : Young Modulus
 %         nu       : Poisson ratio
 %         nodes    : List of nodes
 %         elem     : List of elements of the mesh
 %         order    : order of the elements
 %         gorp     : if 0 : strain at Gauss points
 %                    if 1 : strain at nodes
 %         varargin : if 1 : plane deformations else plane constrains
 
 % output : sigma : stress tensor at nodes or Gauss points

 % First build the model's stiffness matrix
 
 if mat(1) == 0
    E = mat(2); nu = mat(3);
    S = 1/E*[1,-nu,-nu,0,0,0 ; -nu,1,-nu,0,0,0 ; -nu,-nu,1,0,0,0 ;...
             0,0,0,2*(1+nu),0,0 ; 0,0,0,0,2*(1+nu),0 ; 0,0,0,0,0,2*(1+nu)];
    Sm1 = inv(S);
 else
    error('The required material behaviour is not implemented')
 end
 
 if gorp == 0
    sigma = zeros(6*size(elem,1),1);
 else
    sigma = zeros(6*size(nodes,1),1);
 end
 
 for i=1:size(elem,1)
    Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
    nnodes = size(Xloc1,1);
    Xloc = zeros(3*nnodes,1);
    for j=1:nnodes
        mapu(1,[3*j-2,3*j-1,3*j]) = [3*elem(i,j)-2, 3*elem(i,j)-1, 3*elem(i,j)];
        Xloc([3*j-2,3*j-1,3*j],1) = [Xloc1(j,1);Xloc1(j,2);Xloc1(j,3)];
    end
    
    % Build Ke
    Ke = stifmat3D(Xloc,order,Sm1,1);
    % Extract ue
%    mapu = [3*elem(i,1)-2,3*elem(i,1)-1,3*elem(i,1),...
%            3*elem(i,2)-2,3*elem(i,2)-1,3*elem(i,2),...
%            3*elem(i,3)-2,3*elem(i,3)-1,3*elem(i,3),...
%            3*elem(i,4)-2,3*elem(i,4)-1,3*elem(i,4),];
    ue = u(mapu,1);
    
    % Build sigma
    sigmae1 = Ke*ue;
    if gorp == 0
        sigmae  = sigmae1;
        maps = [6*i-5,6*i-4,6*i-3,6*i-2,6*i-1,6*i];
    else                             
        ns = size(sigmae1,1); ne = size(elem,2);
        sigmae = zeros(ns*ne,1); maps = zeros(1,ns*ne);
        for j=1:ne
           sigmae( [1+ns*(j-1):ns*j],1 ) =...
                                sigmae1/ntoelem(elem(i,j),1); % Pass to nodes
           maps(1, [1+ns*(j-1):ns*j]) =...
                              [ 6*elem(i,j)-5,6*elem(i,j)-4,6*elem(i,j)-3, ...
                                6*elem(i,j)-2,6*elem(i,j)-1,6*elem(i,j) ];
        end
    end
    
    sigma(maps,1) = sigma(maps,1) + sigmae;
 end
end