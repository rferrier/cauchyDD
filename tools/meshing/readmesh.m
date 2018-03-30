function [ nodes,elements,ntoelem,boundary,order,physical ] = readmesh( adress )
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
        nodes = fscanf(file,'%*d %f %f %*f',[2 nb])'; % read & store the nodes

    elseif strcmpi(line,'$Elements') % Elements flag
        nb = sscanf(fgetl(file), '%d',1); % number of elements
        elements = zeros(nb,6);
        boundary = zeros(nb,4);
        physical = zeros(nb,1);
        nelem = 1;
        nbound = 1;
        for i=1:nb
            data = sscanf(fgetl(file), '%d');
            if data(2) == 2  % Core elements
                order = 1;   % Yes, that's not very optimized
                elements(nelem,1:3) = data(data(3)+4:end); % Only read index of nodes
                physical(nelem) = data(data(3)+3); % Physical set
                nelem = nelem+1;
            elseif data(2) == 1 % Boundary elements
                boundary(nbound,1:3) = data([data(3)+2,data(3)+4:end]);
                nbound = nbound+1;
            elseif data(2) == 9  % Core elements T6
                order = 2;
                elements(nelem,:) = data(data(3)+4:end); % Only read index of nodes
                physical(nelem) = data(data(3)+3); % Physical set
                nelem = nelem+1;
            elseif data(2) == 8 % Boundary elements S3
                boundary(nbound,:) = data([data(3)+2,data(3)+4:end]);
                nbound = nbound+1;
            else
                error('The dimension of elements you are trying to use is not implemented');
            end
            
        end
        % Cut the database
        if nbound > 1
            elements(nelem:end,:) = [];
            physical(nelem:end) = [];
            boundary(nbound:end,:) = [];
        else
            boundary = [];
        end
        if order == 1
           elements(:,4:end) = [];
           boundary(:,4:end) = [];
        end
    end
 end

 % ntoelem is a table storing the number of elements to witch each node belongs
 ntoelem = zeros(size(nodes,1), 1);
 for i=1:size(elements,1)
     for j=1:size(elements,2)
         ntoelem(elements(i,j),1) = ntoelem(elements(i,j),1) + 1;
     end
 end

end
