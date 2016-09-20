function [index] = findCorner (xl, yl)
 % This function finds the corner of an L-curve

 if size(xl,1) <= 2
    error('Not enough points on L-curve : it has no corner')
 end
 
 % put into loglog
 x = log(xl); y = log(yl);
 
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
 for i=1:sx
    tt = zeros(n+1,1);
    for j=1:n+1
       tt(j) = t(i)^(n+1-j);
    end
    
    % Compute values on this point
    xx(i) = px'*tt; xx1 = px1'*tt; xx2 = px2'*tt;
    yy(i) = py'*tt; yy1 = py1'*tt; yy2 = py2'*tt;
    
    denom = yy2*xx1 - xx2*yy1;
    
    if denom ~= 0
       Rc = (xx1^2 + yy1^2)^(3/2) / denom;
    else
       Rc = R+1.;  % eliminate this candidate
    end
    
    % Test if the candidate is smaller than the current R
    if R == 0 || Rc < R
       R = Rc;
       index = i;
    end
 end
 
end
