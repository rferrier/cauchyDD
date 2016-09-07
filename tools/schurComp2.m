function [ S, b, map ] = schurComp2( K, f, map, mzeros, varargin )
 % This function computes the Schur complement and the Schur Rhs :
 
 ud = zeros(size(K,1),1);
 if numel(varargin)>0
     ud = cell2mat(varargin(1));
 end
 
 % K = [Kii,Kib;Kbi,Kbb]
 % f = [fi;fb]
 
 %map = nindex;%zeros( 1, 2*size(nindex,1) );
 %for i=1:size(nindex,1)
  %   map(1,2*i-1) = 2*nindex(i,1)-1;
   %  map(1,2*i) = 2*nindex(i,1);
 %end
 
 % Assemble Rhs with Dirichlet
 f = f + K*ud;
 
 Kbb = K(map,map);
 fb  = f(map,1);
 
 Kii = K;
 fi  = f;
 Kii([map;mzeros],:) = [];
 Kii(:,[map;mzeros]) = [];
 fi([map;mzeros],:)  = [];
 
 Kbi = K(map,:);
 Kbi(:,[map;mzeros]) = [];
 
 Kib = K(:,map);
 Kib([map;mzeros],:) = [];
 
 S = full(Kbb - Kbi*inv(Kii)*Kib);
 b = fb - Kbi*inv(Kii)*fi;
end