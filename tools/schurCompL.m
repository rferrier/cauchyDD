function [ S, b, map, Iii ] = schurCompL( K, f, nindex, nbloq, c2nodes )
 % This function computes the Schur complement and the Schur Rhs of a
 % matrix with Lagrange multipliers :
 
 % K = [Kii,Kib,Cb';...
 %      Kbi,Kbb,Ci';...
 %      Cb,Ci,0];
 % f = [fi;fb;0]
 
 map = zeros( 1, 2*size(nindex,1) );
 for i=1:size(nindex,1)
     map(1,2*i-1) = 2*nindex(i,1)-1;
     map(1,2*i) = 2*nindex(i,1);
 end
 
 Kbb = K(map,map);
 fb  = f(map,1);
 %Cb  = K(map,c2nodes(map));
 
 Kii = K;
 fi  = f;
 Kii(map,:) = [];
 Kii(:,map) = [];
 fi(map,:)  = [];
 %fi(end-nbloq:end,:) = [];
 
 Kbi = K(map,:);
 Kbi(:,map) = [];
 %Kbi(:,size(Kbi,2)-nbloq:end) = [];
 
 Kib = K(:,map);
 Kib(map,:) = [];
 %Kib(size(Kib,1)-nbloq:end,:) = [];
 
 %nz = 0;  % Debug stuff
 %for i=1:size(Kii,1)
 %    if norm(Kii(:,i))==0
 %        nz = nz+1;
 %    end
 %end
% nz
 
 Iii  = inv(Kii);

 %Iiii = Iii(1:end-nbloq-1,1:end-nbloq-1);

 S = full(Kbb - Kbi*Iii*Kib);
 b = fb - Kbi*Iii*fi;
end