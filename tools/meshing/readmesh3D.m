function [ nodes,elements,ntoelem,boundary,order ] = readmesh3D( adress )
 % Extracts the tables of nodes and elements
 % Adapted from a codebase from Pierre-Eric Allier
 
 file = fopen(adress,'r');
 
 if file <= 0
     error('mesh file not found')
 end
 
 while ~feof(file)  % read until end
    line = fgetl(file);
    if strcmpi(line,'$Nodes') % Nodes flag
        nb = sscanf(fgetl(file), '%d',1); % number of nodes
        nodes = fscanf(file,'%*d %f %f %f',[3 nb])'; % read & store the nodes

    elseif strcmpi(line,'$Elements') % Elements flag
        nb = sscanf(fgetl(file), '%d',1); % number of elements
        elements = zeros(nb,10);
        boundary = zeros(nb,7);
        nelem = 1;
        nbound = 1;
        for i=1:nb
            data = sscanf(fgetl(file), '%d');
            if data(2) == 4  % Core elements
                order = 1;
                elements(nelem,1:4) = data(data(3)+4:end); % Only read index of nodes
                nelem = nelem+1;
            elseif data(2) == 2 % Boundary elements
                boundary(nbound,1:4) = data([data(3)+2,data(3)+4:end]);
                nbound = nbound+1;
            elseif data(2) == 11 % Tetra10 core elements
               order = 2;
               elements(nelem,:) = data(data(3)+4:end); % Only read index of nodes
               nelem = nelem+1;
            elseif data(2) == 9 % Boundary elements (T6)
                boundary(nbound,:) = data([data(3)+2,data(3)+4:end]);
                nbound = nbound+1;
            else
                error('The dimension of elements you are trying to use is not implemented');
            end
            
        end
        % Cut the database
        if nbound > 1
            elements(nelem:end,:) = [];
            boundary(nbound:end,:) = [];
        else
            boundary = [];
        end
        if order == 1
           elements(:,5:end) = [];
           boundary(:,5:end) = [];
        end
    end
 end

 % ntoelem is a table storing the number of elements each node belongs to
 ntoelem = zeros(size(nodes,1), 1);
 for i=1:size(elements,1)
     for j=1:size(elements,2)
         ntoelem(elements(i,j),1) = ntoelem(elements(i,j),1) + 1;
     end
 end

end