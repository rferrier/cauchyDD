function En = energy( u,nodes,elem,mat,order,varargin )
 % This function computes the stresses in the material
 % input : u        : global solution in displacement
 %         E        : Young Modulus
 %         nu       : Poisson ratio
 %         nodes    : List of nodes
 %         elem     : List of elements of the mesh
 %         order    : order of the elements
 %         varargin : if 1 : plane deformations else plane constrains
 
 % mat = [0, E, nu] if isotropic
 % mat = [1, Ex, Ey, nuxy, Gxy] TODO : add rotation of the axis
 % mat = [1.5, k11, k12, k22, k33] orthotropic, but with identification coefficients

 % output : sigma : stress tensor at nodes or Gauss points

 % First build the model's stiffness matrix
 mo = 0;
 if numel(varargin)>0
     mo = cell2mat(varargin(1));
 end
 
 % Test if the material is homogeneous
 if size(elem,2)~=3 && size(elem,2)~=6 % TODO : more stable system
    inh = 1;
    physical = elem(:,1);
    elem = elem(:,2:end);
    if min(size(mat))==0
       warning('First material is empty');
    end
    nmat = size(mat,2);
 else
    inh = 0;
    physical = ones(size(elem,1));
    mat = {mat};
    nmat = 1;
 end

 Sm1 = {};
 for i=1:nmat
    matt = mat{i};
    if matt(1) == 0
       E = matt(2); nu = matt(3);
       if mo == 0
          S = 1/E*[1,-nu,0 ; -nu,1,0 ; 0,0,2*(1+nu)];
          Sm1{i} = inv(S);
       else
          kappa = E/(3*(1-2*nu));
          mu = E/(2*(1+nu));
          Sm1{i} = [4*mu/3+kappa, kappa-2*mu/3, 0;...
                 kappa-2*mu/3, 4*mu/3+kappa, 0;...
                 0, 0, mu];
       end
    elseif matt(1) == 1
       Ex = matt(2); Ey = matt(3); nuxy = matt(4); Gxy = matt(5);
       S = [1/Ex, -nuxy/Ex, 0;...  % Rem : nuxy/Ex = nuyx/Ey
            -nuxy/Ex, 1/Ey, 0;...
            0, 0, 1/Gxy];
       Sm1{i} = inv(S);
       if mo == 1   % TODO : difference between plan constraints and deformations.
          warning('Using plane deformation, not plane constraints as required')
       end
    elseif matt(1) == 1.5
       k11 = matt(2); k12 = matt(3); k22 = matt(4); k33 = matt(5);
       Sm1{i} = [k11, k12, 0;... 
            k12, k22, 0;...
            0, 0, k33];
       if mo == 1   % TODO : difference between plan constraints and deformations.
          warning('Using plane deformation, not plane constraints as required')
       end
    else
       error('The required material behaviour is not implemented')
    end
 end
 
 En = zeros(size(nodes,1),1);

 for i=1:size(elem,1)
    phys = physical(i);
    Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
    nnodes = size(Xloc1,1);
    Xloc = zeros(2*nnodes,1);
     for j=1:nnodes
         mapu(1,[2*j-1,2*j]) = [2*elem(i,j)-1, 2*elem(i,j)];
         Xloc([2*j-1,2*j],1) = [Xloc1(j,1);Xloc1(j,2)];
     end
    
    % Build Ke
    Ke = stifmat(Xloc,order,Sm1{phys},0);

    % Extract ue
    ue = u(mapu,1);

    % Build e
    En(i) = ue'*Ke*ue;
 end
 
end
