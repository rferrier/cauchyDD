function [indm,p] = findPicard3 (x, varargin)
 % This function finds the points in a solution vector, that don't respect the 
 % discrete Picard condition.
 % This function interpolates the vector x with a polynomial
 % and finds the point from witch the derivate is > 0

 %warning('off','all','local'); % Because of high order singularity issue

 npo = size(x,1);
 if npo <= 2
    error(['Not enough points on curve : it has no corner', ...
      ' Rem : if you actually have lots of points, try to transpose the vectors'])
 end
 
 d = 0;   % relation between nb of dof and size of the vector
 if numel(varargin)>0
     d = cell2mat(varargin(1));
 end

 t    = 1:npo;     % Dummy axis
 n    = floor(npo/d)+1;  % degree of polynoms
 %n    = min(n,20); % Not too high degree

 p = polyfit(t,x',n)';
 
 tt = zeros(n+1,npo);
 for j=1:n+1
    tt(j,:) = t.^(n+1-j);
 end

 px = p'*tt;  % interpolation

 inferior = [];

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
 
 % build the derivation matrix
 Md = zeros(n+1);
 for i=2:n+1
    Md(i,i-1) = n-i+2;
 end
    
 p1 = Md*p;  % Derivative
 px1 = p1'*tt; % interpolation of the derivative
   
 posit = find( px1>0 ); % Find the first point with px1>0
 if size(posit,2)>0
    indm1 = posit(1);
 end

 pxp = px; pxp(inferior) = max(px)+1; % Makes sure inferior values won't be selected
 noposit = setdiff(1:size(px,2),posit);
 pxp(noposit) = max(px)+1; % Makes sure negative derivative won't be selected
    
 if size(noposit,2)>0 % If there is a negative part
    nopo1 = noposit(1);
    pxp(1:nopo1) = max(px)+1; % Makes sure first values won't be selected
 end

 [minpx,indm1] = min(pxp);
     
 if size(posit,2)==0 % There is no minimum : select the last one
    indm = npo;
 else
    % Now, we got our minimum, find the maximum before that minimum
    [maxpx,indm2] = max(px(1:indm1));
    % And find the first value over that maximum
    cand = find( x(indm1+1:end)>= maxpx);
    indm = indm1+min(cand)-1;
    if min(size(indm))==0
       indm = npo;
    end
 end

end
