function epsilon = strain( u,nodes,elem,order,gorp,ntoelem )
 % This function computes the stresses in the material
 % input : u        : global solution in displacement
 %         nodes    : List of nodes
 %         elem     : List of elements of the mesh
 %         order    : order of the elements
 %         gorp     : if 0 : strain at Gauss points
 %                    if 1 : strain at nodes
 
 % output : epsilon : stress tensor at nodes or Gauss points
 
 if gorp == 0
    epsilon = zeros(3*size(elem,1),1);
 else
    epsilon = zeros(3*size(nodes,1),1);
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
    Ke = strainmat(Xloc,order,Sm1,1);

    % Extract ue
    ue = u(mapu,1);
    
    % Build epsilon
    epsilone1 = Ke*ue;
    if gorp == 0
        epsilone  = epsilone1;
        maps = [3*i-2,3*i-1,3*i];
    else
        ns = size(epsilone1,1); ne = size(elem,2);
        epsilone = zeros(ns*ne,1); maps = zeros(1,ns*ne);
        for j=1:ne
           epsilone( [1+ns*(j-1):ns*j],1 ) =...
                                epsilone1/ntoelem(elem(i,j),1); % Pass to nodes
           maps(1, [1+ns*(j-1):ns*j]) =...
                              [3*elem(i,j)-2,3*elem(i,j)-1,3*elem(i,j)];
        end
    end
    
    epsilon(maps,1) = epsilon(maps,1) + epsilone;
 end
 
end
