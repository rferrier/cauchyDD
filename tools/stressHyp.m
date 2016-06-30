function sigma = stressHyp( u,E,nu,nodes,elem,order,gorp,ntoelem,alpha,varargin )
 % This function computes the stresses in the material
 % input : u        : global solution in displacement
 %         E        : Young Modulus
 %         nu       : Poisson ratio
 %         nodes    : List of nodes
 %         elem     : List of elements of the mesh
 %         order    : order of the elements
 %         gorp     : if 0 : strain at Gauss points
 %                    if 1 : strain at nodes
 %                    if 2 : compute residual
 %         ntoelem  : nb of elements each node is linked to
 %         alpha    : compressibility parameter
 %         varargin : if 1 : plane deformations else plane constrains
 
 % output : sigma : stress tensor at nodes or Gauss points

 % First build the model's stiffness matrix
 mo = 0;
 if numel(varargin)>0
     mo = cell2mat(varargin(1));
 end
 
 if gorp == 0
    sigma = zeros(3*size(elem,1),1);
 elseif gorp == 1
    sigma = zeros(3*size(nodes,1),1);
 elseif gorp == 2
    sigma = zeros(2*size(nodes,1),1);
 end
 
 for i=1:size(elem,1)
    Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
    nnodes = size(Xloc1,1);
    Xloc = zeros(2*nnodes,1);
    for j=1:nnodes
        Xloc([2*j-1,2*j],1) = [Xloc1(j,1);Xloc1(j,2)];
    end
    
    if mo == 0
        S = 1/E*[1,-nu,0 ; -nu,1,0 ; 0,0,1+nu];
        Sm1 = inv(S);
    else
        kappa = E/(3*(1-2*nu));
        mu = E/(2*(1+nu));

        Sm1 = [4*mu/3+kappa, kappa-2*mu/3, 0;...
               kappa-2*mu/3, 4*mu/3+kappa, 0;...
               0, 0, 2*mu];
    end
    
    % Build Ke
    Ke = stifmat(Xloc,order,Sm1,1);
    % Extract ue
    mapu = [2*elem(i,1)-1,2*elem(i,1),2*elem(i,2)-1,...
           2*elem(i,2),2*elem(i,3)-1,2*elem(i,3)];
    ue = u(mapu,1);

    % Compute non-linea term
    nnodes = size(Xloc1,1);
    uloc = zeros(2*nnodes,1);      % Extract local displacement
    for j=1:nnodes
        uloc([2*j-1,2*j],1) = [ u( 2*elem(i,j)-1 ) ; u( 2*elem(i,j) ) ];
    end

    x11 = Xloc(1,1); x21 = Xloc(3,1); x31 = Xloc(5,1);  % nodal coordinates
    x12 = Xloc(2,1); x22 = Xloc(4,1); x32 = Xloc(6,1);
    S = abs( .5*((x21-x11)*(x32-x12)-(x31-x11)*(x22-x12)) );% element area
    Be=[x22-x32,0,x32-x12,0,x12-x22,0;...
        0,x31-x21,0,x11-x31,0,x21-x11;
        x31-x21,x22-x32,x11-x31, ...
        x32-x12,x21-x11,x12-x22]/(2*S);
    % Compute tr(epsilon)
     epsilon = Be*uloc;
     treps = epsilon(1)+epsilon(2);
    
    % Build sigma
    sigmae1 = Ke*ue + alpha*kappa*treps^3*[1;1;0];
    
    if gorp == 0
        sigmae  = sigmae1;
        maps = [3*i-2,3*i-1,3*i];
        
    elseif gorp == 1
        
        sigmae = [sigmae1/ntoelem(elem(i,1),1);
                  sigmae1/ntoelem(elem(i,2),1);
                  sigmae1/ntoelem(elem(i,3),1)]; % pass to Gauss points
              
        maps = [3*elem(i,1)-2,3*elem(i,1)-1,3*elem(i,1),...
                3*elem(i,2)-2,3*elem(i,2)-1,3*elem(i,2),...
                3*elem(i,3)-2,3*elem(i,3)-1,3*elem(i,3)];
            
    elseif gorp == 2
        
        x11 = Xloc(1,1); x21 = Xloc(3,1); x31 = Xloc(5,1);  % nodal coordinates
        x12 = Xloc(2,1); x22 = Xloc(4,1); x32 = Xloc(6,1);
        S = abs( .5*((x21-x11)*(x32-x12)-(x31-x11)*(x22-x12)) );% element area
        Be=[x22-x32,0,x32-x12,0,x12-x22,0;...
            0,x31-x21,0,x11-x31,0,x21-x11;
            x31-x21,x22-x32,x11-x31, ...
            x32-x12,x21-x11,x12-x22]/(2*S);
        
        % This is not exactly the stress, but the RHS
        sigmae  = S*Be'*sigmae1;
        maps = [2*elem(i,1)-1,2*elem(i,1),...
                2*elem(i,2)-1,2*elem(i,2),...
                2*elem(i,3)-1,2*elem(i,3)];
    end
    sigma(maps,1) = sigma(maps,1) + sigmae;
 end
 
end

