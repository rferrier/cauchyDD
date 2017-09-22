function ret = plotGMSH3D( fields, elem, nodes, name )
 % This function writes an output file readable with Gmsh
 % Adaptated from a codebase from ?? (someone I don't know)
 
% prepares output file
fmid = fopen(['meshes/',name,'.msh'],'w');
fprintf(fmid,'%s\n%s\n%s\n','$MeshFormat','2.2 0 8','$EndMeshFormat');
nnodes = size(nodes,1);

for i=1:size(fields,1)
    field = cell2mat(fields(i,1));
    nb = size(field,1);
    legend = cell2mat(fields(i,2));
    % First : what is the kind of field ?

    for j=1:size(field,2)
        if nb == nnodes
            % scalar output
            fprintf(fmid,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n','$NodeData','1',['"',legend,'"'],'1','0.0','3',num2str(j-1),'1');
            fprintf(fmid,'%d\n',nnodes);
            for n=1:nnodes
                fprintf(fmid,'%d %E \n',n,field(n,j));   
            end
            fprintf(fmid,'%s\n','$EndNodeData');

        elseif nb == 3*nnodes
            % vector output
            fprintf(fmid,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n','$NodeData','1',['"',legend,'"'],'1','0.0','3',num2str(j-1),'3');
            fprintf(fmid,'%d\n',nnodes);
            for n=1:nnodes
                fprintf(fmid,'%d %E %E %E\n',n,field(3*n-2,j),field(3*n-1,j),field(3*n,j));   
            end
            fprintf(fmid,'%s\n','$EndNodeData');

        elseif nb == 6*nnodes
            % tensor output
            fprintf(fmid,'%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n','$NodeData','1',['"',legend,'"'],'1','0.0','3',num2str(j-1),'9');
            fprintf(fmid,'%d\n',nnodes);
            for n=1:nnodes
                fprintf(fmid,'%d %E %E %E %E %E %E %E %E %E\n',n,...
                      field(6*n-5,j),field(6*n-2,j),field(6*n-1,j),field(6*n-2,j),field(6*n-4,j),...
                      field(6*n,j),field(6*n-1,j),field(6*n,j),field(6*n-3,j));   
            end
            fprintf(fmid,'%s\n','$EndNodeData');

        else
            error(' Error at output : could not find the type of field')
        end
    end
end

fclose(fmid);

end