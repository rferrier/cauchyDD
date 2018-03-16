function dist = angleDist( angle1, angle2 )
% This function computes a distance between 2 vectors of angles,
% taking into account the fact that the angles are periodic.

 % Use the dots corresponding to the angles
 dot1 = [ cos(angle1), sin(angle1) ];
 dot2 = [ cos(angle2), sin(angle2) ];

 dist = (dot1(:,1)-dot2(:,1)).^2 + (dot1(:,2)-dot2(:,2)).^2;
 dist = norm(dist);
end
