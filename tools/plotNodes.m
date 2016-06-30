function [ret, retdefo] = plotNodes( u,elems,nodes )
 % This function plots a field defined at the nodes of the mesh
 % Adaptated from a codebase from Pierre-Eric Allier

 % Plot mesh
 ret = patch('Faces',elems(:,1:3),'Vertices',nodes,'FaceAlpha',0);

 nodesdefo = nodes + 1e2*reshape(u,2,[])';
 retdefo = patch('Faces',elems(:,1:3),'Vertices',nodesdefo,'FaceAlpha',0);
 
 % TODO : use the code on r>2013
%  set(ret,'FaceAlpha',1);
% 
%  nnodes = size(nodes,1);
%  c = uicontextmenu;
%  N = numel(u)/nnodes;
%  if N > 1
%      for i=1:N
%          uimenu(c,'Label',['u' num2str(i)],'Callback',{@setdata,ret,i});
%      end
%  end
%  u1 = reshape(u,N,[])';
%  set(ret,'UserData',u1);
%  set(ret,'FaceVertexCData',u1(:,1));
%  set(ret,'FaceColor','interp');
%  set(ret,'LineStyle','none');
%  % On Maltab r>2013, the leatest line shouldn't return a warning
%  set(ret,'UIContextMenu',c);

end
