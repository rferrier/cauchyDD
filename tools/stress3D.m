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
             0,0,0,1+nu,0,0 ; 0,0,0,0,1+nu,0 ; 0,0,0,0,0,1+nu];
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
        Xloc([3*j-2,3*j-1,3*j],1) = [Xloc1(j,1);Xloc1(j,2);Xloc1(j,3)];
    end
    
    % Build Ke
    Ke = stifmat3D(Xloc,order,Sm1,1);
    % Extract ue
    mapu = [3*elem(i,1)-2,3*elem(i,1)-1,3*elem(i,1),...
            3*elem(i,2)-2,3*elem(i,2)-1,3*elem(i,2),...
            3*elem(i,3)-2,3*elem(i,3)-1,3*elem(i,3),...
            3*elem(i,4)-2,3*elem(i,4)-1,3*elem(i,4),];
    ue = u(mapu,1);
    
    % Build sigma
    sigmae1 = Ke*ue;
    if gorp == 0
        sigmae  = sigmae1;
        maps = [6*i-5,6*i-4,6*i-3,6*i-2,6*i-1,6*i];
    else
        sigmae = [sigmae1/ntoelem(elem(i,1),1);
                  sigmae1/ntoelem(elem(i,2),1);
                  sigmae1/ntoelem(elem(i,3),1);
                  sigmae1/ntoelem(elem(i,4),1)]; % pass to Gauss points
              
        maps = [6*elem(i,1)-5,6*elem(i,1)-4,6*elem(i,1)-3,...
                            6*elem(i,1)-2,6*elem(i,1)-1,6*elem(i,1),...
                6*elem(i,2)-5,6*elem(i,2)-4,6*elem(i,2)-3,...
                            6*elem(i,2)-2,6*elem(i,2)-1,6*elem(i,2),...
                6*elem(i,3)-5,6*elem(i,3)-4,6*elem(i,3)-3,...
                             6*elem(i,3)-2,6*elem(i,3)-1,6*elem(i,3),...
                6*elem(i,4)-5,6*elem(i,4)-4,6*elem(i,4)-3,...
                             6*elem(i,4)-2,6*elem(i,4)-1,6*elem(i,4),];
    end
    
    sigma(maps,1) = sigma(maps,1) + sigmae;
 end
 
end

