function [indm,p] = findPicard2 (x, varargin)
 % This function finds the points in a solution vector, that don't respect the 
 % discrete Picard condition.
 % This function interpolates the vector x with a polynomial
 % and finds the point from witch the derivate is > 0
 npo = size(x,1);
 if npo <= 2
    error(['Not enough points on curve : it has no corner', ...
      ' Rem : if you actually have lots of points, try to transpose the vectors'])
 end
 
 d = 0;   % relation between nb of dof and size of the vector
 if numel(varargin)>0
     d = cell2mat(varargin(1));
 end
 
 maxi = 0;   % detect only the max
 if numel(varargin)>1
     maxi = cell2mat(varargin(2));
 end
 
 crit = 1;   % Criterion : 1=derivative, 2=minimum
 if numel(varargin)>2
     crit = cell2mat(varargin(3));
 end
 
 t    = 1:npo;     % Dummy axis
 n  = floor(npo/d)+1;  % degree of polynoms
 indm = 0;

 p = polyfit(t,x',n)';
 
 tt = zeros(n+1,npo);
 for j=1:n+1
    tt(j,:) = t.^(n+1-j);
 end

 px = p'*tt;  % interpolation

 inferior = [];
 if maxi == 1  % suppress the inferior values
    inferior = find( x<px' );
    tprim = t; tprim(inferior) = [];
    xprim = x; xprim(inferior) = [];
    n  = min( max(2,floor(size(xprim,1)/d+1)), size(xprim,1) );  % reactualize degree of polynoms
    p = polyfit(tprim,xprim',n)';
    
    tt = zeros(n+1,npo);
    for j=1:n+1
       tt(j,:) = t.^(n+1-j);
    end

    px = p'*tt;  % interpolation
 end
 
 
 if crit == 1
    % build the derivation matrix
    Md = zeros(n+1);
    for i=2:n+1
       Md(i,i-1) = n-i+2;
    end
    
    p1 = Md*p;  % Derivative
    px1 = p1'*tt; % interpolation of the derivative
   
    posit = find( px1>0 ); % Find the first point with px1>0
    if size(posit>0)
       indm = posit(1);
    end
 elseif crit == 2
    px(inferior) = max(px)+1; % Makes sure inverior values von't be selected
    [px,indm] = min(px);
%    indm = find ( px == min(px) );
 else
    error('Third optional argument is not readable')
 end
 
end