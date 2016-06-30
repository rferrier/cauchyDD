function M = H1demi( nb, nodes, boundary, entity )
%H1demi computes the H1/2 mass matrix on a boundary

 % First term
 inodes = zeros(size(nodes,1), 1);
 h = 0; % Total length of the boundary
 N1 = zeros(nb);
 for i=1:size(boundary)
     if boundary(i,1) == entity
        if inodes(boundary(i,2)) == 0
            no1 = boundary(i,2);
            inodes(no1) = 1; %Don't re-compute this one
            inodes2 = zeros(size(nodes,1), 1);
            inodes2(no1) = 1; % sum i \neq j
            for j=1:size(boundary)
                if boundary(j,1) == entity
                    if inodes2(boundary(j,2)) == 0
                        no2 = boundary(j,2);
                        x1 = nodes(no1,1); y1 = nodes(no1,2);
                        x2 = nodes(no2,1); y2 = nodes(no2,2);
                        map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
                        len = sqrt( (x2-x1)^2 + (y2-y1)^2 );
                        h = max(h, len);  % Re-actualise h

                        inodes2(no2) = 1;
                        N1(map,map) = N1(map,map) + ...
                            1/(len^2)*[ 1,0,-1,0; 0,1,0,-1;...
                                        -1,0,1,0; 0,-1,0,1];
                    end
                    if inodes2(boundary(j,3)) == 0
                        no2 = boundary(j,3);
                        x1 = nodes(no1,1); y1 = nodes(no1,2);
                        x2 = nodes(no2,1); y2 = nodes(no2,2);
                        map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
                        len = sqrt( (x2-x1)^2 + (y2-y1)^2 );
                        h = max(h, len);  % Re-actualise h

                        inodes2(no2) = 1;
                        N1(map,map) = N1(map,map) + ...
                            1/(len^2)*[ 1,0,-1,0; 0,1,0,-1;...
                                        -1,0,1,0; 0,-1,0,1];
                    end
                end
            end
        end

        % Second node of the element
        if inodes(boundary(i,3)) == 0
            no1 = boundary(i,3);
            inodes(no1) = 1;
            inodes2 = zeros(size(nodes,1), 1);
            inodes2(no1) = 1;
            for j=1:size(boundary)
                if boundary(j,1) == entity
                    if inodes2(boundary(j,2)) == 0
                        no2 = boundary(j,2);
                        x1 = nodes(no1,1); y1 = nodes(no1,2);
                        x2 = nodes(no2,1); y2 = nodes(no2,2);
                        map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
                        len = sqrt( (x2-x1)^2 + (y2-y1)^2 );
                        h = max(h, len);  % Re-actualise h

                        inodes2(no2) = 1;
                        N1(map,map) = N1(map,map) + ...
                            1/(len^2)*[ 1,0,-1,0; 0,1,0,-1;...
                                        -1,0,1,0; 0,-1,0,1];
                    end
                    if inodes2(boundary(j,3)) == 0
                        no2 = boundary(j,3);
                        x1 = nodes(no1,1); y1 = nodes(no1,2);
                        x2 = nodes(no2,1); y2 = nodes(no2,2);
                        map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
                        len = sqrt( (x2-x1)^2 + (y2-y1)^2 );
                        h = max(h, len);  % Re-actualise h

                        inodes2(no2) = 1;
                        N1(map,map) = N1(map,map) + ...
                            1/(len^2)*[ 1,0,-1,0; 0,1,0,-1;...
                                        -1,0,1,0; 0,-1,0,1];
                    end
                end
            end
        end
     end
 end

 % Second term
 N2 = zeros(nb);
 for i=1:size(boundary)
     if boundary(i,1) == entity
         no1 = boundary(i,2); no2 = boundary(i,3);
         x1 = nodes(no1,1); y1 = nodes(no1,2);
         x2 = nodes(no2,1); y2 = nodes(no2,2);
         len = sqrt((x1-x2)^2 + (y1-y2)^2);
         map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
         N2(map,map) = N2(map,map) + len/2*[ 1,0,1/3,0; 0,1,0,1/3;...
                                  1/3,0,1,0; 0,1/3,0,1];
     end
 end

 M = h*N1 + N2/h^2;

end

