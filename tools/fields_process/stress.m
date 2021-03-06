function sigma = stress( u,E,nu,nodes,elem,order,gorp,ntoelem,varargin )
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
 mo = 0;
 if numel(varargin)>0
     mo = cell2mat(varargin(1));
 end
 
 if mo == 0
    S = 1/E*[1,-nu,0 ; -nu,1,0 ; 0,0,2*(1+nu)];
    Sm1 = inv(S);
 else
    kappa = E/(3*(1-2*nu));
    mu = E/(2*(1+nu));
    Sm1 = [4*mu/3+kappa, kappa-2*mu/3, 0;...
           kappa-2*mu/3, 4*mu/3+kappa, 0;...
           0, 0, mu];
 end
 
 if gorp == 0
    sigma = zeros(3*size(elem,1),1);
 else
    sigma = zeros(3*size(nodes,1),1);
 end
 
 for i=1:size(elem,1)
    Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
    nnodes = size(Xloc1,1);
    Xloc = zeros(2*nnodes,1);
     for j=1:nnodes
         mapu(1,[2*j-1,2*j]) = [2*elem(i,j)-1, 2*elem(i,j)];
         Xloc([2*j-1,2*j],1) = [Xloc1(j,1);Xloc1(j,2)];
     end
    
    % Build Ke
    Ke = stifmat(Xloc,order,Sm1,1);

    % Extract ue
    ue = u(mapu,1);
    
    % Build sigma
    sigmae1 = Ke*ue;
    if gorp == 0
        sigmae  = sigmae1;
        maps = [3*i-2,3*i-1,3*i];
    else
        ns = size(sigmae1,1); ne = size(elem,2);
        sigmae = zeros(ns*ne,1); maps = zeros(1,ns*ne);
        for j=1:ne
           sigmae( [1+ns*(j-1):ns*j],1 ) =...
                                sigmae1/ntoelem(elem(i,j),1); % Pass to nodes
           maps(1, [1+ns*(j-1):ns*j]) =...
                              [3*elem(i,j)-2,3*elem(i,j)-1,3*elem(i,j)];
        end
    end
    
    sigma(maps,1) = sigma(maps,1) + sigmae;
 end
 
end
