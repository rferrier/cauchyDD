% 15/11/2017 : Correction d'un champ par minimisation de la dérivée
% Version modernisée
 
clear all;
close all;
 
Norm       = 1;
ordpD      = 15;
ordp       = 15;
sideL      = 20; % Nb of points on the domain (for integral approximation)
pm4        = 4; % Position of the 1d line
upper_term = 0; % If 0 : keep only terms st i+j<ordp
id_crack   = 0;  % Should I identify the crack (binary stuff) ?
threshold  = .1;  % Threshold for the crack

%VAR = load('fields/AnBnoplanel615.mat');
%VAR = load('fields/AnBsmile615.mat');
%VAR = load('fields/AnBlosange2l615.mat');
%VAR = load('fields/AnBlosange615.mat');
%VAR = load('fields/AnBm15b2.mat');
%VAR = load('fields/AnBs1013.mat');
VAR = load('fields/AnBm15.mat');
%VAR = load('fields/AnBs15.mat');
%VAR = load('fields/AnBs15b1.mat');
%VAR = load('fields/AnBs215b1.mat');
%VAR = load('fields/AnBs215b1.mat');
Lhso = VAR.Lhs; Rhso = VAR.Rhs; Xs = VAR.Xs; Ys = VAR.Ys; ordpo = VAR.ordp;
X0 = VAR.X0; Y0 = VAR.Y0; Lx = VAR.Lx; Ly = VAR.Ly;
L1x = VAR.L1x; L2x = VAR.L2x; L1y = VAR.L1y; L2y = VAR.L2y;
uplo1 = VAR.uplo1; uplo = VAR.uplo;

% Re-build the smaller Lhs and Rhs
Lhs = zeros((ordp+1)^2); Rhs = zeros((ordp+1)^2,1);
for k = 0:ordp
   for l = 0:ordp
      Rhs(1+l+(ordp+1)*k) = Rhso(1+l+(ordpo+1)*k);
      for i = 0:ordp
         for j = 0:ordp
            Lhs( 1+l+(ordp+1)*k , 1+j+(ordp+1)*i ) = ...
                              Lhso( 1+l+(ordpo+1)*k , 1+j+(ordpo+1)*i );
         end
      end
   end
end

if upper_term == 0 % Remove upper order terms
   toremove = [];
   for i=0:ordp
      for j=0:ordp
         if i+j>ordp
            toremove(end+1) = 1+j+(ordp+1)*i;
         end
      end
   end
   tokeep = setdiff(1:size(Rhs,1),toremove);
   
   toremoveD = [];
   for i=0:ordpD
      for j=0:ordpD
         if i+j>ordpD
            toremoveD(end+1) = 1+j+(ordpD+1)*i;
         end
      end
   end
   tokeepD = setdiff(1:(ordpD+1)^2,toremoveD);
else
   tokeep  = 1:size(Rhs,1);
   tokeepD = 1:(ordpD+1)^2;
end

nxs = (max(Xs)-min(Xs))/100; nys = (max(Ys)-min(Ys))/100;
X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
% The reference
figure;
hold on;
surf(X,Y,uplo1);
shading interp;
colorbar();
axis('equal');

McCoef = zeros(size(Rhs));
McCoef(tokeep) = -Lhs(tokeep,tokeep)\Rhs(tokeep);
solup = zeros(101,101);
for k=0:ordp
   for l=0:ordp
      solup = solup + McCoef(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
   end
end
solup = solup';
% Find the coefficients on the rest of the base, that minimize the derivative.
% First, build a basis that is orthogonal (in sense of Legendre) to the first vectors

% Matrix for the L2 Inner product
PS0 = zeros((ordpD+1)^2);
for i=0:ordpD
   for j=0:ordpD
      for k=0:ordpD
         for l=0:ordpD
            ordx = i+k+1;
            ordy = j+l+1;
            %if mod(ordx,1) == 0 && mod(ordy,1) == 0 % Provided L1x = -L2x && L1y = -L2y
               PS0(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) = ...
                    Lx^2*(L2x^ordx - L1x^ordx)/ordx * (L2y^ordy - L1y^ordy)/ordy;
                    % Lx* beacuse there is a variable change x' = Lx*x and y'=Lx*y
            %end
         end
      end
   end
end

%figure
%plot(log(abs(eig(PS0))));

% Reshape the matrix
DiagMult = sqrt(diag(1./diag(PS0)));

%PS0p = zeros((ordpD+1)^2); % Re-compute with the shaped form (actually, I'm not sure it's useful)
%for i=0:ordpD
%   for j=0:ordpD
%      for k=0:ordpD
%         for l=0:ordpD
%            ordx = i+k+1;
%            ordy = j+l+1;
%            if mod(ordx,1) == 0 && mod(ordy,1) == 0 % Provided L1x = -L2x && L1y = -L2y
%               PS0p(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) = ...
%                    DiagMult(j+1+(ordpD+1)*i,j+1+(ordpD+1)*i) * ...
%                    DiagMult(l+1+(ordpD+1)*k,l+1+(ordpD+1)*k) * ...
%                    Lx^2*(L2x^ordx - L1x^ordx)/ordx * (L2y^ordy - L1y^ordy)/ordy;
%                    % Lx* beacuse there is a variable change x' = Lx*x and y'=Lx*y
%            end
%         end
%      end
%   end
%end

PS0p = DiagMult'*PS0*DiagMult; % (This reshape appears to be useless)
%PS0p = PS0;

%    [PassD, ~] = eig(LhsO); % This matrix gives the Legendre basis
%    PassD = eye((ordpD+1)^2);
[PS0p,dist] = closestSPD(PS0p); % Just in case PS0p is not SPD because of floats
if dist > 1e-5
   warning('Frobenius distance to closest SPD is',num2str(dist));
end
PassD1 = chol(PS0p,'lower');     %CHOLESKI
%for i=1:size(PassD1,2) % Normalize the stuff
%   PassD1(:,i) = PassD1(:,i) / norm(PassD1(:,i));
%end
PassD = inv(PassD1');
PassD = DiagMult * PassD ;

%PS0D = PassD'*PS0*PassD;  % Debug check
%bug;

tic
 
if Norm == 2
   % Matrix for the "H1-L2" Inner product
   PS1 = zeros((ordpD+1)^2);
   for i=0:ordpD
      for j=0:ordpD
         for k=0:ordpD
            for l=0:ordpD
               if k+i>1
                  PS1(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) = ...
                        PS1(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) + ...
                        1/Lx^2*k*i*(L2x^(k+i-1)-L1x^(k+i-1))*(L2y^(l+j+1)-L1y^(l+j+1))/...
                                             ((k+i-1)*(l+j+1));
               end
               if l+j>1
                  PS1(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) = ...
                        PS1(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) +...
                        1/Ly^2*l*j*(L2x^(k+i+1)-L1x^(k+i+1))*(L2y^(l+j-1)-L1y^(l+j-1))/...
                                             ((k+i+1)*(l+j-1));
               end
            end
         end
      end
   end
   
   RhsR = zeros((ordpD+1)^2,1);
   for i=0:ordp
      for j=0:ordp
         for k=0:ordpD
            for l=0:ordpD
               if i+k>1
                  RhsR(l+1+(ordpD+1)*k) = RhsR(l+1+(ordpD+1)*k) - McCoef(1+j+(ordp+1)*i) * ...
                        1/Lx^2*i*k*(L2x^(i+k-1)-L1x^(i+k-1))*(L2y^(j+l+1)-L1y^(j+l+1))/((i+k-1)*(j+l+1));
               end
               if j+l>1
                  RhsR(l+1+(ordpD+1)*k) = RhsR(l+1+(ordpD+1)*k) - McCoef(1+j+(ordp+1)*i) * ...
                        1/Ly^2*j*l*(L2x^(i+k+1)-L1x^(i+k+1))*(L2y^(j+l-1)-L1y^(j+l-1))/((i+k+1)*(j+l-1));
               end
            end
         end
      end
   end

   RhsD = zeros((ordpD+1)^2,1); % Will be truncated later : for now, it contains all the coeffs from 0 to ordpD

%   if upper_term == 0 % Put upper terms to 0
%      RhsR(toremoveD) = 0;
%      PS1(toremoveD,:) = 0;
%      PS1(:,toremoveD) = 0;
%%      PS1(toremoveD,toremoveD) = eye(size(toremove,2));
%   end
   
   LhsD = PassD'*PS1*PassD;
   RhsD1 = PassD'*RhsR;
   % We are only using the last vectors from the basis : put zeros
   for k=0:ordp
      for l=0:ordp
%          if (k<=ordp && l<=ordp)
            LhsD(l+1+(ordpD+1)*k,:) = 0;
            LhsD(:,l+1+(ordpD+1)*k) = 0;
%         end
      end
   end

   % At this point, there are plenty of zero columns (and lines) in LhsD.
   % We are going to add identity on those in order to make the system invertible
   % This has no impact on the solution as the RhsD in front is 0.
   for i=0:ordpD
      for j=0:ordpD
         if norm(LhsD(j+1+(ordpD+1)*i,:)) == 0
            LhsD(j+1+(ordpD+1)*i,j+1+(ordpD+1)*i) = 1;
            RhsD1(j+1+(ordpD+1)*i) = 0;
         end
      end
   end
   RhsD = RhsD1;
   McCoefD = LhsD\RhsD;
   
   % Go back to the regular basis (not the orthogonal one)
   McCoefD0 = PassD*McCoefD;

   % Add it to solup
   X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
   solupp = solup'; % Because of the transpose
%   for k=0:ordpD
%      for l=0:ordpD
%         for ii=0:k
%            for jj=0:l
%               solupp = solupp + PassD(jj+1+(ordpD+1)*ii, l+1+(ordpD+1)*k) * ...
%                       McCoefD(1+l+(ordpD+1)*k) .* (X/Lx-X0)'.^ii * (Y/Lx-Y0).^jj;
%            end
%         end
%      end
%   end
   for k=0:ordpD
      for l=0:ordpD
         solupp = solupp + ...
                  McCoefD0(1+l+(ordpD+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
      end
   end
   solupp = solupp';

   figure;
   hold on;
   surf(X,Y,solupp);
   shading interp;
   colorbar();
   axis('equal');
   
   error = norm( solupp(:)-uplo1(:) ) / norm(uplo1(:));
   
   if id_crack == 1 % Identify the crack
      crackID = zeros(size(solupp));
      maxto = max(max(solupp))*threshold;
      tofind = find( solupp>maxto )-1;
      indicei = rem( tofind , size(solupp,1) )+1;
      indicej = floor( tofind / size(solupp,1) )+1;
      
      for i = 1:size(tofind,1)
         crackID(indicei(i),indicej(i)) = 1;
      end
      
      crackRef = zeros(size(uplo1));
      maxto = max(max(uplo1))*1e-5;
      tofind = find( uplo1>maxto )-1;
      indicei = rem( tofind , size(uplo1,1) )+1;
      indicej = floor( tofind / size(uplo1,1) )+1;
      
      for i = 1:size(tofind,1)
         crackRef(indicei(i),indicej(i)) = 1;
      end
   
      figure;
      hold on;
      surf(X,Y,crackRef);
      shading interp;
      colorbar();
      axis('equal');
   
      figure;
      hold on;
      surf(X,Y,crackID);
      shading interp;
      colorbar();
      axis('equal');
      
      errorID = norm(crackID(:)-crackRef(:))/norm(crackRef(:));
   end

elseif Norm == 1 % L1 norm of the gradient (new beautyful stuff)
   % First task : build the matrix giving the derivatives from the coefficients
   nxs = (max(Xs)-min(Xs))/sideL; nys = (max(Ys)-min(Ys))/sideL;
   X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
    
   Ax0 = zeros((sideL+1)^2,(ordpD+1)^2); Ay0 = zeros((sideL+1)^2,(ordpD+1)^2);
   for k=1:ordpD
      for l=0:ordpD
         Amute = 1/Lx*k*(X/Lx-X0)'.^(k-1) * (Y/Lx-Y0).^l;
         Ax0(:,l+1+(ordpD+1)*k) = Amute(:);
         % Lx* beacuse there is a variable change x' = Lx*x and y'=Lx*y
      end
   end
   for k=0:ordpD
      for l=1:ordpD
         Amute = 1/Ly*(X/Lx-X0)'.^k * l*(Y/Lx-Y0).^(l-1);
         Ay0(:,l+1+(ordpD+1)*k) = Amute(:);
         % Lx* beacuse there is a variable change x' = Lx*x and y'=Lx*y
      end
   end
    
   Ax = Ax0*PassD; % Put it in the "Legendre" basis
   Ay = Ay0*PassD;
    
   % ... And re-write over them for the first coefficients (those that are not in the Legendre basis)
   for k=1:ordp
      for l=0:ordp
         Amute = 1/Lx*k*(X/Lx-X0)'.^(k-1) * (Y/Lx-Y0).^l;
         Ax(:,l+1+(ordpD+1)*k) = Amute(:);
      end
   end
   for k=0:ordp
      for l=1:ordp
         Amute = 1/Ly*(X/Lx-X0)'.^k * l*(Y/Lx-Y0).^(l-1);
         Ay(:,l+1+(ordpD+1)*k) = Amute(:);
      end
   end
   
   % Build Isum with the 1/2 on the borders and 1/4 on the corners
   Is2 = ones(sideL+1,1); Is2(1) = .5; Is2(end) = .5;
   Isum = Is2*Is2'; Isum = Isum(:);
   
   % Build the matrix of the linear optim problem
   C = [ zeros(1,(ordpD+1)^2) , Isum' , Isum' ]; %ones(1,2*(sideL+1)^2)

   % Now, the constraint matrix
   Im = zeros(0,(ordpD+1)^2);
   Id = eye((sideL+1)^2);
   Zoo = zeros((sideL+1)^2); Zno = zeros( (sideL+1)^2 , (ordpD+1)^2 );
   
%   Itor(toremove,toremove) = eye(size(toremove)); % For the 
   ind = 1;
   for k=0:ordp
      for l=0:ordp
         Im( ind , 1+l+(ordpD+1)*k ) = 1; % Already-determined coefs
         ind = ind+1;
      end
   end

   Zoom = zeros(size(Im),(sideL+1)^2);
   
   Ac = [ Im, Zoom, Zoom ;...
          Ax, -Id, Zoo ;...
          -Ax, -Id, Zoo ;...
          Ay, Zoo, -Id ;...
          -Ay, Zoo, -Id];
          
   % Lhs of the constraint
   bc = zeros( 4*(sideL+1)^2 , 1 );
   
   a0 = zeros( 0 , 1 ); % Already determined coefficients
   ind = 1;
   for k=0:ordp
      for l=0:ordp
         a0(ind,1) = McCoef(l+1+(ordp+1)*k,1); 
         ind = ind+1;
      end
   end
   
   bc = [a0;bc];
   
   lb = -1e5*ones(2*(sideL+1)^2+(ordpD+1)^2,1);
   ub = 1e5*ones(2*(sideL+1)^2+(ordpD+1)^2,1); % Bounds
   
   cvect = 'C'*ones( 2*(sideL+1)^2 + (ordpD+1)^2 , 1 ); cvect = char(cvect);
   ctype1 = 'S'*ones( ind-1 , 1 );
   ctype2 = 'U'*ones( 4*(sideL+1)^2 , 1 );
   ctype = char( [ctype1;ctype2] );
   
%   if upper_term == 0
%      warning('upper terms not implemented for L1 minimisation');
%   end
   
   %% Call the Solver
   param.lpsolver = 1;
   [ sol2 , minim , errnum , extra ] = glpk ( C' , Ac , bc , lb , ub , ctype , cvect , 1 , param); %(1:27,:)
   sol2 = sol2(1:(ordpD+1)^2);
   
%   %% Alternative : homemade stuff
%   alp = Im\a0;
%   N = null(Im);
%   if size(N,2) > 0
%      sol0 = zeros(size(N,2),1);
%      fctest0 = Isum'*abs(Ax*(alp+N*sol0)) + Isum'*abs(Ay*(alp+N*sol0));
%      ix0 = Ax*(alp+N*sol0); iy0 = Ay*(alp+N*sol0);
%      
%      sol = sol0;
%      solp = sol;
%      
%      Sxr = zeros( size(Ax*N*sol), 100); Syr = zeros( size(Ay*N*sol), 100);
%      Xr  = zeros( size(Ax*N*sol), 100); Yr  = zeros( size(Ay*N*sol), 100);
%      
%      relax = 1;
%      
%      i = 1;
%      for j=1:1000
%         % Compute signs
%         ix = Ax*(alp+N*solp); iy = Ay*(alp+N*solp);
%         Sx = sign(ix); Sy = sign(iy);
%         Xr(:,i) = ix; Yr(:,i) = iy;
%         Sxr(:,i) = Sx; Syr(:,i) = Sy;
%         
%         % Right hand side
%         b = transpose( (Isum.*Sx)'*Ax*N + (Isum.*Sy)'*Ay*N );
%         b = transpose( (Isum.*Sx.*abs(ix0))'*Ax*N + (Isum.*Sy.*abs(iy0))'*Ay*N );
%         b = transpose( (Isum.*ix)'*Ax*N + (Isum.*iy)'*Ax*N );
%         b = transpose( (Isum.*(ix/mean(abs(ix))))'*Ax*N + (Isum.*(iy/mean(abs(iy))))'*Ax*N );
%         
%         dire = b; % Steepest descent direction
%         
%         % Randomize
%         bcand = zeros( size(N,2) , 10*size(N,2)+1 ); socand = zeros(size(bcand));
%         fcand = zeros(1,10*size(N,2)+1);
%         for k=1:10*size(N,2)
%            bcand(:,k)  = abs( rand(size(N,2),1) ) .* dire;
%            socand(:,k) = solp - relax*bcand(:,k);
%            fcand(k) = Isum'*abs(Ax*(alp+N*socand(:,k))) + Isum'*abs(Ay*(alp+N*socand(:,k)));
%         end
%         bcand(:,end) = dire;
%         socand(:,end) = solp - relax*bcand(:,end);
%         fcand(end) = Isum'*abs(Ax*(alp+N*socand(:,end))) + Isum'*abs(Ay*(alp+N*socand(:,end)));
%         
%         % Choose the best one
%         [~,icand] = min(fcand);
%         dire = bcand(:,icand);
%         
%         sol = solp - relax*dire;
%         relaxrec(i)  = relax;
%         fctest(i)    = Isum'*abs(Ax*(alp+N*sol)) + Isum'*abs(Ay*(alp+N*sol));
%         
%         if i>1
%            difftest(i) = fctest(i)-fctest(i-1);
%         else
%            difftest(i) = fctest(i) - fctest0; % First is 0
%         end
%         
%         if i>1 
%            crit2(i) = norm(solp-sol)/norm(solp);
%            ratio(i) = difftest(i)/norm(solp-sol);
%         else
%            crit2(i) = 1;
%         end
%         if crit2(i) < 1e-10
%            break;
%         end
%         
%         indexfirst = max(1,i-1);
%         criterion = max(difftest(indexfirst:i));
%         if criterion <= 0 && relax <= .5
%            relax = 2*relax;
%            relaxmax = relax;
%         elseif difftest(i) > 0
%            relax = .5*relax;
%            i = i-1;
%            sol = solp; % Cancel the iteration
%         end
%            
%         solp = sol;
%         i = i+1;
%      end
%      
%      sol = alp + N*sol;
%   else
%      sol = alp;
%   end
%   
%   bug;
%   %%
%   %% Alternative : homemade stuff
%   alp = Im\a0;
%   N = null(Im);
%   if size(N,2) > 0
%      sol0 = zeros(size(N,2),1); 
%%      X0 = [Ax*N,zeros(size(Ax,1),size(N,2)) ; ...
%%            zeros(size(Ay,1),size(N,2)), Ay*N]*[sol0;sol0];
%      X0 = [Ax*N ; Ay*N]*sol0;
%
%      fctest0 = Isum'*abs(Ax*(alp+N*sol0)) + Isum'*abs(Ay*(alp+N*sol0));
%      ix0 = Ax*(alp+N*sol0); iy0 = Ay*(alp+N*sol0);
%      
%%      sol = sol0; X  = [Ax*N,zeros(size(Ax,1),size(N,2)) ; ...
%%                        zeros(size(Ay,1),size(N,2)), Ay*N]*[sol;sol];
%%      solp = sol; Xp = [Ax*N,zeros(size(Ax,1),size(N,2)) ; ...
%%                        zeros(size(Ay,1),size(N,2)), Ay*N]*[solp;solp];
%      sol = sol0; X  = [Ax*N;Ay*N]*sol;
%      solp = sol; Xp = [Ax*N;Ay*N]*solp;
%      
%      Sxr = zeros( size(Ax*N*sol), 100); Syr = zeros( size(Ay*N*sol), 100);
%      Xr  = zeros( size(Ax*N*sol), 100); Yr  = zeros( size(Ay*N*sol), 100);
%      
%      relax = 1;
%      
%      i = 1;
%      for j=1:100
%         % Compute signs
%         ix = Ax*(alp+N*solp); iy = Ay*(alp+N*solp);
%         Sx = sign(ix); Sy = sign(iy);
%         Xr(:,i) = ix; Yr(:,i) = iy;
%         Sxr(:,i) = Sx; Syr(:,i) = Sy;
%         
%         % Right hand side
%%         b = transpose( (Isum.*Sx)'*Ax*N + (Isum.*Sy)'*Ay*N );
%         b = transpose( [(Isum.*Sx)', (Isum.*Sy)'] );
%%         b = transpose( (Isum.*Sx.*abs(ix0))'*Ax*N + (Isum.*Sy.*abs(iy0))'*Ay*N );
%%         b = transpose( (Isum.*ix)'*Ax*N + (Isum.*iy)'*Ax*N );
%%         b = transpose( (Isum.*(ix/mean(abs(ix))))'*Ax*N + (Isum.*(iy/mean(abs(iy))))'*Ax*N );
%         
%         dire = b; % Steepest descent direction
%         
%%         % Randomize
%%         bcand = zeros( size(N,2) , 10*size(N,2)+1 ); socand = zeros(size(bcand));
%%         fcand = zeros(1,10*size(N,2)+1);
%%         for k=1:10*size(N,2)
%%            bcand(:,k)  = abs( rand(size(N,2),1) ) .* dire;
%%            socand(:,k) = solp - relax*bcand(:,k);
%%            fcand(k) = Isum'*abs(Ax*(alp+N*socand(:,k))) + Isum'*abs(Ay*(alp+N*socand(:,k)));
%%         end
%%         bcand(:,end) = dire;
%%         socand(:,end) = solp - relax*bcand(:,end);
%%         fcand(end) = Isum'*abs(Ax*(alp+N*socand(:,end))) + Isum'*abs(Ay*(alp+N*socand(:,end)));
%%         
%%         % Choose the best one
%%         [~,icand] = min(fcand);
%%         dire = bcand(:,icand);
%         
%         X = Xp - relax*dire;
%         relaxrec(i)  = relax;
%         sol = [Ax*N; Ay*N]\X;  % Projection on the right field
%         fctest(i)    = Isum'*abs(Ax*(alp+N*sol)) + Isum'*abs(Ay*(alp+N*sol));
%         
%         if i>1
%            difftest(i) = fctest(i)-fctest(i-1);
%         else
%            difftest(i) = fctest(i) - fctest0; % First is 0
%         end
%         
%         if i>1 
%            crit2(i) = norm(solp-sol)/norm(solp);
%            ratio(i) = difftest(i)/norm(solp-sol);
%         else
%            crit2(i) = 1;
%         end
%         if crit2(i) < 1e-10
%            break;
%         end
%         
%         indexfirst = max(1,i-1);
%         criterion = max(difftest(indexfirst:i));
%         if criterion <= 0 && relax <= .5
%            relax = 2*relax;
%            relaxmax = relax;
%         elseif difftest(i) > 0
%            relax = .5*relax;
%            i = i-1;
%            sol = solp; % Cancel the iteration
%            X = Xp;
%         end
%            
%         solp = sol; Xp = [Ax*N;Ay*N]*solp;
%         i = i+1;
%      end
%      
%      sol = alp + N*sol;
%   else
%      sol = alp;
%   end
%   
%   bug;
   %%

   sol = sol2;
   
   % Recover the result
%   McCoefDt = zeros(4*(ordpD+1)^2,1);
%   McCoefDt(~toremove) = sol;
%   McCoefDt(~toremove) = sol;

   McCoefDt  = sol(1:(ordpD+1)^2);
%   Xx        = sol((ordpD+1)^2+1:(ordpD+1)^2+(sideL+1)^2);
%   Yy        = sol((ordpD+1)^2+(sideL+1)^2+1:end);
   
   % Remove the already existing coeffs (as it will be added later)
   McCoefD = McCoefDt;
   McCoef0 = McCoefDt-McCoefDt;
   for k=0:ordp
      for l=0:ordp
         McCoefD(l+1+(ordpD+1)*k,1) = 0;
         McCoef0(l+1+(ordpD+1)*k,1) = McCoef(l+1+(ordp+1)*k,1);
      end
   end
   
%   % Debug : see the value without any McCoefDt
%   Xx0 = abs(Ax*McCoef0); Yy0 = abs(Ay*McCoef0);
%   val0 = C*[McCoef0;Xx0;Yy0];
   
   % Add it to solup
   nxs = (max(Xs)-min(Xs))/100; nys = (max(Ys)-min(Ys))/100;
   X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
   solupp = solup'; % Because of the transpose
   for k=0:ordpD
      for l=0:ordpD
         for ii=0:k
            for jj=0:l
               solupp = solupp + PassD(jj+1+(ordpD+1)*ii, l+1+(ordpD+1)*k) * ...
                       McCoefD(1+l+(ordpD+1)*k) .* (X/Lx-X0)'.^ii * (Y/Lx-Y0).^jj;
            end
         end
      end
   end
   solupp = solupp';
   
   figure;
   hold on;
   surf(X,Y,solupp);
   shading interp;
   colorbar();
   axis('equal');
   
   error = norm( solupp(:)-uplo1(:) ) / norm(uplo1(:));
   
   if id_crack == 1 % Identify the crack
      crackID = zeros(size(solupp));
      maxto = max(max(solupp))*threshold;
      tofind = find( solupp>maxto )-1;
      indicei = rem( tofind , size(solupp,1) )+1;
      indicej = floor( tofind / size(solupp,1) )+1;
      
      for i = 1:size(tofind,1)
         crackID(indicei(i),indicej(i)) = 1;
      end
      
      crackRef = zeros(size(uplo1));
      maxto = max(max(uplo1))*1e-5;
      tofind = find( uplo1>maxto )-1;
      indicei = rem( tofind , size(uplo1,1) )+1;
      indicej = floor( tofind / size(uplo1,1) )+1;
      
      for i = 1:size(tofind,1)
         crackRef(indicei(i),indicej(i)) = 1;
      end
   
      figure;
      hold on;
      surf(X,Y,crackRef);
      shading interp;
      colorbar();
      axis('equal');
   
      figure;
      hold on;
      surf(X,Y,crackID);
      shading interp;
      colorbar();
      axis('equal');
      
      errorID = norm(crackID(:)-crackRef(:))/norm(crackRef(:));
   end
   
else % No gradient stuff
   solupp = solup;
   McCoefD = zeros((ordpD+1)^2,1);
end
 
toc
 
% Plot on the line X = 4
figure;
hold on;
nys = (max(Ys)-min(Ys))/100;
Y = min(Ys):nys:max(Ys); X = pm4;

solup = 0;
for k=0:ordp
   for l=0:ordp
      solup = solup + McCoef(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
   end
end

%plot( Y, solup, 'Color', 'black' );

solupp = solup;
for k=0:ordpD
   for l=0:ordpD

      for ii=0:k
         for jj=0:l
            solupp = solupp + PassD(jj+1+(ordpD+1)*ii, l+1+(ordpD+1)*k) * ...
                    McCoefD(1+l+(ordpD+1)*k) .* (X/Lx-X0)'.^ii * (Y/Lx-Y0).^jj;
         end
      end
  end
end
plot( Y, solupp, 'Color', 'green' );

% Reference
plot( Y, uplo(3:3:end), 'Color', 'red' );
