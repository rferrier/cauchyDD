function mesh2GMSH( nodes, elements, boundary, name )
% This function exports a mesh to gmsh
     
 nnodes = size(nodes,1); nelem = size(elements,1); nbound = size(boundary,1);

 fmid = fopen(['meshes/',name,'.msh'],'w');
 fprintf(fmid,'%s\n%s\n%s\n','$MeshFormat','2.2 0 8','$EndMeshFormat');
     
 % Nodes
 fprintf(fmid,'%s\n','$Nodes');
 fprintf(fmid,'%d\n',nnodes);

 if size(nodes,2) == 2
    for n=1:nnodes
        fprintf(fmid,'%d %d %d %d \n',n,nodes(n,1),nodes(n,2),0);   
    end
 else % There is a third one
    for n=1:nnodes
        fprintf(fmid,'%d %d %d %d \n',n,nodes(n,1),nodes(n,2),nodes(n,3));   
    end
 end
 fprintf(fmid,'%s\n','$EndNodes');
    
 % Elements
 fprintf(fmid,'%s\n','$Elements');
 fprintf(fmid,'%d\n',nelem+nbound);
 for n=1:nbound
     fprintf(fmid,'%d %d %d %d %d %d %d \n',n,1,2,boundary(n,1),...
         boundary(n,1),boundary(n,2),boundary(n,3));   
 end
 for n=1:nelem
     fprintf(fmid,'%d %d %d %d %d %d %d %d \n',n+nbound,2,2,5,1,...
         elements(n,1),elements(n,2),elements(n,3));   
 end
 fprintf(fmid,'%s\n','$EndElements');
 
 fclose(fmid);
end
