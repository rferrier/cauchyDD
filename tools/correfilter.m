function [b] = correfilter (nodes, Lc, method, varargin)
% This function does a truncation of a signal with the Lc correlation length.
% nodes    : nodes of the mesh to apply the white noise toascii
% Lc       : correlation length
% method   : 'rectangle'
% varargin : stochastic germ (vector) or signal to filter

nnodes = size(nodes, 1);

 % Manage the user defined random vector
 if numel(varargin)>0
    germ = cell2mat(varargin(1));
 else
    germ = rand( 2*nnodes, 1 );
 end

 % Build the Correlation matrix
 Mc  = zeros(2*nnodes);
 nno = zeros(nnodes,1); % Stores the number of nodes
 if method == 'rectangle'  % /!\ Matrix won't be symmetric
    for i=1:nnodes
       xi = nodes(i,1);  yi = nodes(i,2);
       for j=1:nnodes
          xj = nodes(j,1);  yj = nodes(j,2);
          dist2 = (xi-xj)^2 + (yi-yj)^2;
          if dist2 <= Lc^2
             Mc(2*i-1,2*j-1) = 1;
             Mc(2*i,2*j) = 1;
             nno(i) = nno(i)+1;
          end
       end
       % normalize
       Mc(2*i-1,:) = Mc(2*i-1,:)/nno(i);
       Mc(2*i,:) = Mc(2*i,:)/nno(i);
    end
    
 else % Unknown method
    error(['The method ', method, ' for white noise generation is unimplemented']);
 end

 b = Mc*germ;
end
