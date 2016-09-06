function [ S, b, map ] = schurComp2( K, f, map, zeros )
 % This function computes the Schur complement and the Schur Rhs :
 
 % K = [Kii,Kib;Kbi,Kbb]
 % f = [fi;fb]
 
 %map = nindex;%zeros( 1, 2*size(nindex,1) );
 %for i=1:size(nindex,1)
  %   map(1,2*i-1) = 2*nindex(i,1)-1;
   %  map(1,2*i) = 2*nindex(i,1);
 %end
 
 Kbb = K(map,map);
 fb  = f(map,1);
 
 Kii = K;
 fi  = f;
 Kii([map;zeros],:) = [];
 Kii(:,[map;zeros]) = [];
 fi([map;zeros],:)  = [];
 
 Kbi = K(map,:);
 Kbi(:,[map;zeros]) = [];
 
 Kib = K(:,map);
 Kib([map;zeros],:) = [];
 
 S = full(Kbb - Kbi*inv(Kii)*Kib);
 b = fb - Kbi*inv(Kii)*fi;
end