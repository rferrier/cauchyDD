function E = elasticEnergy (u,mat,nodes,elem,order,varargin)
% This function returns the elastic energy corresponding to a displacement
% field in a domain

 mo = 0;
 if numel(varargin)>0
     mo = cell2mat(varargin(1));
 end
 
 if mat(1) == 0
    E = mat(2); nu = mat(3);
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
 elseif mat(1) == 1
    Ex = mat(2); Ey = mat(3); nuxy = mat(4); Gxy = mat(5);
    S = [1/Ex, -nuxy/E1, 0;...
         -nuxy/E1, 1/Ey, 0;...
         0, 0, 1/Gxy];
    if mo == 1   % TODO : difference between plan constraints and deformations.
       warning('Using plane deformation, not plane constraints as required')
    end
 elseif mat(1) == 2
    % Random inhomogeneous medium
    E = mat(2); nu = mat(3); minE = mat(4); maxE = mat(5);
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
 else
    error('The required material behaviour is not implemented')
 end
 
 E = 0;
 
 for i=1:size(elem,1)
    Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
    nnodes = size(Xloc1,1);
    Xloc = zeros(2*nnodes,1);
     for j=1:nnodes
         mapu(1,[2*j-1,2*j]) = [2*elem(i,j)-1, 2*elem(i,j)];
         Xloc([2*j-1,2*j],1) = [Xloc1(j,1);Xloc1(j,2)];
     end
    
    % Build Ke
    Ke = stifmat(Xloc,order,Sm1,0);

    % Extract ue
    ue = u(mapu,1);
    
    % Increment E
    E = E + ue'*Ke*ue;
 end

end