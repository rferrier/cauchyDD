function [indm,a1,b1,a2,b2] = findPicard (x, varargin)
 % This function finds the points in a solution vector, that don't respect the 
 % discrete Picard condition.
 % This function tries to interpolate the vector x with 2 lines going from
 % [0;t_1] and the other from [t_1;t_end] : f1=a1x+b1 and f2 = a2x+b2
 % the t_1 minimizing the total MS error is the right one
 npo = size(x,1);
 if npo <= 2
    error(['Not enough points on curve : it has no corner', ...
      ' Rem : if you actually have lots of points, try to transpose the vectors'])
 end
 
 remove = 0;   % number of nodes you have to remove
 if numel(varargin)>0
     remove = cell2mat(varargin(1));
 end
 
 t    = 1:npo;     % Dummy axis
 emin = 2*norm(x); % Initialize with something greater than any error
 indm = 0;
 a1 = 0; b1 = 0; a2 = 0; b2 = 0;
 
 for i=2+remove:npo-1-remove % try all the possible corners
    p1   = polyfit(t(1:i),x(1:i)',1);  px1 = p1(2) + p1(1)*t;
    err1 = norm(px1(1:i)-x(1:i)');
    p2 = polyfit(t(i:end),x(i:end)',1);  % i and not i+1
    px2 = p2(2) + p2(1)*t; err2 = norm(px2(i:end)-x(i:end)');
    erto = err1+err2; % /!\ point i counts 2 times
    
    if erto <= emin % store this one
       emin = erto; indm = i;
       a1 = p1(1); b1 = p1(2); a2 = p2(1); b2 = p2(2);
    end
    
%    hold on;
%    plot(x,'Color','black')
%    plot(px1-x,'Color','yellow')
%    plot(px2-x,'Color','green')
%    legend('solution coefficients','left interpolation', 'right interpolation')
%    figure;
    
 end

end