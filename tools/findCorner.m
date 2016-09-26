function [index] = findCorner (xl, yl)
 % This function finds the corner of an L-curve

 if size(xl,1) <= 2
    error('Not enough points on L-curve : it has no corner')
 end
 
 % put into loglog
 x = log10(xl); y = log10(yl);
 
 sx = size(x,1);  % size of x (hopefully too the size of y)
 n  = floor(sx/2)+1;%max(sx-2,floor(sx/2)+1);   % degree of polynoms
 t  = (1:1:sx)';   % coarse mesh
 tp = (1:.1:sx)';  % fine mesh (useless for now)
 
 % First, interpolate x and y
 px = polyfit(t,x,n)';
 py = polyfit(t,y,n)';
 
 % build the derivation matrix
 Md = zeros(n+1);
 for i=2:n+1
    Md(i,i-1) = n-i+2;
 end
 
 % Then derivate
 px1 = Md*px;
 px2 = Md*px1;
 py1 = Md*py;
 py2 = Md*py1;
 
 % Find the point with the smallest curve radius
 R = 0;
 index = 0;
 xx = zeros(sx,1);  yy = zeros(sx,1); %debug stuff
 xx1 = zeros(sx,1);  yy1 = zeros(sx,1);
 xx2 = zeros(sx,1);  yy2 = zeros(sx,1);
 Ga = zeros(sx,1);
 for i=2:sx-1  % First and last are forbitten (because of interpolation bound effects)
    tt = zeros(n+1,1);
    for j=1:n+1
       tt(j) = t(i)^(n+1-j);
    end
    
    % Compute values on this point
    xx(i) = px'*tt; xx1(i) = px1'*tt; xx2(i) = px2'*tt;
    yy(i) = py'*tt; yy1(i) = py1'*tt; yy2(i) = py2'*tt;
    
    denom = (xx1(i)^2 + yy1(i)^2)^(3/2);
    
    if denom ~= 0
       Ga(i) =  (yy2(i)*xx1(i) - xx2(i)*yy1(i)) / denom;
    else
       Ga(i) = R+1.;  % eliminate this candidate
    end
    
    % Test if the candidate is smaller than the current curvature
    if R == 0 || Ga(i) < R % Ga < 0
       R = Ga(i);
       index = i;
    end
 end
 
end
