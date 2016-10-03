function [b] = correfilter (nodes, Lc, method, varargin)
% This function does a truncation of a signal with the Lc correlation length.
% nodes    : nodes of the mesh to apply the white noise toascii
% Lc       : correlation length
% method   : frequential window : 'rectangle', 'triangle'/'Bartlett',
%                                 'Hann'/'Hanning', 'Hamming',
%                                 'Blackmann'
% varargin : 1- stochastic germ (vector) or signal to filter
%            2- if 0, no division ie this keeps the amplitude of the signal

nnodes = size(nodes, 1);

 % Manage the user defined random vector
 if numel(varargin)>0
    germ = cell2mat(varargin(1));
 else
    germ = randn( 2*nnodes, 1 );
 end
 
 if numel(varargin)>1
    divide = cell2mat(varargin(2));
 else
    divide = 1;
 end

 % Build the Correlation matrix
 Mc  = zeros(2*nnodes);
 nno = zeros(nnodes,1); % Stores the number of nodes
 if strcmp(method, 'rectangle')  % /!\ Matrix won't be symmetric
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
       if divide == 1
          Mc(2*i-1,:) = Mc(2*i-1,:)/nno(i);
          Mc(2*i,:) = Mc(2*i,:)/nno(i);
       end
    end
    
 elseif strcmp(method, 'triangle') || strcmp(method, 'Bartlett')
    for i=1:nnodes
       xi = nodes(i,1);  yi = nodes(i,2);
       for j=1:nnodes
          xj = nodes(j,1);  yj = nodes(j,2);
          dist2 = (xi-xj)^2 + (yi-yj)^2;
          if dist2 <= Lc^2
             toadd = ( Lc - sqrt(dist2) ) / Lc;
             Mc(2*i-1,2*j-1) = toadd;
             Mc(2*i,2*j) = toadd;
             nno(i) = nno(i) + toadd;
          end
       end
       % normalize
       if divide == 1
          Mc(2*i-1,:) = Mc(2*i-1,:)/nno(i);
          Mc(2*i,:) = Mc(2*i,:)/nno(i);
       end
    end
    
 elseif strcmp(method, 'Hann') || strcmp(method, 'Hanning')
    for i=1:nnodes
       xi = nodes(i,1);  yi = nodes(i,2);
       for j=1:nnodes
          xj = nodes(j,1);  yj = nodes(j,2);
          dist2 = (xi-xj)^2 + (yi-yj)^2;
          if dist2 <= Lc^2
             toadd = .5 + .5*cos( pi*sqrt(dist2)/Lc );
             Mc(2*i-1,2*j-1) = toadd;
             Mc(2*i,2*j) = toadd;
             nno(i) = nno(i) + toadd;
          end
       end
       % normalize
       if divide == 1
          Mc(2*i-1,:) = Mc(2*i-1,:)/nno(i);
          Mc(2*i,:) = Mc(2*i,:)/nno(i);
       end
    end
    
 elseif strcmp(method, 'Hamming')
    for i=1:nnodes
       xi = nodes(i,1);  yi = nodes(i,2);
       for j=1:nnodes
          xj = nodes(j,1);  yj = nodes(j,2);
          dist2 = (xi-xj)^2 + (yi-yj)^2;
          if dist2 <= Lc^2
             toadd = .54 + .46*cos( pi*sqrt(dist2)/Lc );
             Mc(2*i-1,2*j-1) = toadd;
             Mc(2*i,2*j) = toadd;
             nno(i) = nno(i) + toadd;
          end
       end
       % normalize
       if divide == 1
          Mc(2*i-1,:) = Mc(2*i-1,:)/nno(i);
          Mc(2*i,:) = Mc(2*i,:)/nno(i);
       end
    end
    
 elseif strcmp(method, 'Blackmann')
    for i=1:nnodes
       xi = nodes(i,1);  yi = nodes(i,2);
       for j=1:nnodes
          xj = nodes(j,1);  yj = nodes(j,2);
          dist2 = (xi-xj)^2 + (yi-yj)^2;
          if dist2 <= Lc^2
             toadd = .42 + .5*cos( pi*sqrt(dist2)/Lc ) ...
                     - .08*cos( 2*pi*sqrt(dist2)/Lc ) ;
             Mc(2*i-1,2*j-1) = toadd;
             Mc(2*i,2*j) = toadd;
             nno(i) = nno(i) + toadd;
          end
       end
       % normalize
       if divide == 1
          Mc(2*i-1,:) = Mc(2*i-1,:)/nno(i);
          Mc(2*i,:) = Mc(2*i,:)/nno(i);
       end
    end
    
 else % Unknown method
    error(['The method ', method,...
            ' for white noise generation is unimplemented']);
 end

 b = Mc*germ;
end
