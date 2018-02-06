% 30/11/2017 : Identification par minimisation d'une fonctionnelle
% régularisée par Total Variation (steepest régularisé)
 
clear all;
close all;

sideL      = 10; % Nb of points on the domain (for integral approximation)
pm4        = -4; % Position of the 1d line
mu         = 1e-4;%1e-4;%5e-4;%1e-3;%.7e-2;%5e-4;%1e-4;%5e-5; % Regularization parameter
%mu         = 1e-1;%1e-2;
epsilon    = 1e-7; % Regularization parameter for TV
ordp       = 9;
normU      = 1;   % Use L1 or L2 minimization
jmax       = 20;  % Coefficient for Picard stuff (if normU == 3)
upper_term = 0;
id_crack   = 0;  % Should I identify the crack (binary stuff) ?
threshold  = .1;  % Threshold for the crack

%VAR = load('fields/AnBnoplanel615.mat');
%VAR = load('fields/AnBs15.mat');
VAR = load('fields/AnBm15b2.mat');
%VAR = load('fields/AnBm15.mat');
%VAR = load('fields/AnBs15.mat');
%VAR = load('fields/AnBs1013.mat');
%VAR = load('fields/AnBs15b1.mat');
%VAR = load('fields/AnBlosangel615.mat');
%VAR = load('fields/AnBlosange2l615.mat');
%VAR = load('fields/AnBsmile615.mat');
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

if upper_term == 0
   toremove = [];
   for i=0:ordp
      for j=0:ordp
         if i+j>ordp
            toremove(end+1) = 1+j+(ordp+1)*i;
         end
      end
   end
   tokeep = setdiff(1:size(Rhs,1),toremove);
%   Lhs(toremove,:) = 0; Lhs(:,toremove) = 0;
%   Rhs(toremove) = 0;
else
   tokeep = 1:size(Rhs,1);
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

% Orthogonality is useless in this case
tic

% First task : build the matrix giving the derivatives from the coefficients
nxs = (max(Xs)-min(Xs))/sideL; nys = (max(Ys)-min(Ys))/sideL;
X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
Ax = zeros((sideL+1)^2,(ordp+1)^2); Ay = zeros((sideL+1)^2,(ordp+1)^2); % Derivative
Bxy = zeros((sideL+1)^2,(ordp+1)^2); % Value

for k=0:ordp
   for l=0:ordp
      Bmute = (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
      Bxy(:,l+1+(ordp+1)*k) = Bmute(:);
   end
end
for k=1:ordp
   for l=0:ordp
      Amute = 1/Lx*k*(X/Lx-X0)'.^(k-1) * (Y/Lx-Y0).^l;
      Ax(:,l+1+(ordp+1)*k) = Amute(:);
   end
end
for k=0:ordp
   for l=1:ordp
      Amute = 1/Ly*(X/Lx-X0)'.^k * l*(Y/Lx-Y0).^(l-1);
      Ay(:,l+1+(ordp+1)*k) = Amute(:);
   end
end
   
% Build Isum with the 1/2 on the borders and 1/4 on the corners
Is2 = ones(sideL+1,1); Is2(1) = .5; Is2(end) = .5;
Isum = Is2*Is2'; Isum = Isum(:);
Isum = Isum*nxs*nys;
    
%%% Newton-Raphson
%if normU == 1
%   sol   = zeros(size(Rhs),1); % Initialize
%   crit  = zeros(100,1); crit2 = zeros(100,1); diff = zeros(100,1);
%   sol0 = -Lhs\Rhs;
%   
%   Sxr = zeros( size(Ax*sol), 100); Syr = zeros( size(Ay*sol), 100);
%   Xr  = zeros( size(Ax*sol), 100); Yr  = zeros( size(Ay*sol), 100); 
%
%   fctest = zeros(100,1); difftest = zeros(100,1);
%   residu = zeros(size(Rhs),100); nres = zeros(100,1);
%   
%   sol = sol0; 
%   solp = sol;
%                  
%   fctest0 = .5*sol'*Lhs*sol + sol'*Rhs +...
%                  mu*(Isum'*sqrt((Ax*sol).^2+epsilon)+Isum'*sqrt((Ay*sol).^2+epsilon));
%   
%   i = 1;
%   for j=1:10
%      % Compute signs
%      ix = Ax*solp; iy = Ay*solp;
%%      Sx = sign(ix); Sy = sign(iy);
%      Sx = ix./(sqrt(ix.^2)+epsilon);   Sy = iy./(sqrt(iy.^2)+epsilon);
%      Ux = epsilon/(ix.^2+epsilon).^(3/2); Uy = epsilon/(iy.^2+epsilon).^(3/2);
%      Xr(:,i) = ix; Yr(:,i) = iy;
%      
%      Sxr(:,i) = Sx; Syr(:,i) = Sy;
%      
%      % Right hand side
%%      b = mu*transpose( (Isum.*Sx)'*Ax + (Isum.*Sy)'*Ay );
%
%%      resid = Lhs*solp + b + Rhs; % Gradient (residual)
%      resid = Lhs*solp + Rhs + mu*( Ax'*( Isum.*Sx) + Ay'*(Isum.*Sy) ); % Gradient (residual)
%      residu(:,i) = resid; nres(i) = norm(resid);
%      Atan = Lhs + mu*( Ax'*(Isum.*Ux)*(Isum.*Ux)'*Ax +...
%                        Ay'*(Isum.*Uy)*(Isum.*Ux)'*Ay );
%
%      sol = solp - Atan\resid;
%      
%      quadratic(i) = .5*sol'*Lhs*sol + sol'*Rhs;
%      absolute(i)  = mu*(Isum'*sqrt((Ax*sol).^2+epsilon)+Isum'*sqrt((Ay*sol).^2+epsilon));
%      fctest(i)    = quadratic(i) + absolute(i);
%
%      solp = sol;
%      
%      if nres(i) < 1e-10 % Stopping condition
%         break;
%      end
%      
%      i = i+1;
%   end
%   
%   figure;
%   plot(log10(nres(1:i)));
   
%% Relaxed Newton %Steepest descent
if normU == 1
   sol   = zeros(size(Rhs),1); % Initialize
   crit  = zeros(100,1); crit2 = zeros(100,1); diff = zeros(100,1);
   relax = 1;
   sol0 = sol;
   sol0(tokeep) = -Lhs(tokeep,tokeep)\Rhs(tokeep);
   
   Sxr = zeros( size(Ax*sol), 100); Syr = zeros( size(Ay*sol), 100);
   Xr  = zeros( size(Ax*sol), 100); Yr  = zeros( size(Ay*sol), 100); 
   relaxrec = zeros(100,1); % Store the relax parameter (debug)
   fctest = zeros(100,1); difftest = zeros(100,1);
   
   sol = sol0; 
   solp = sol;
                  
   fctest0 = .5*sol'*Lhs*sol + sol'*Rhs +...
                  mu*(Isum'*sqrt((Ax*sol).^2+epsilon)+Isum'*sqrt((Ay*sol).^2+epsilon));
   
   i = 1;
   for j=1:5000
      % Compute signs
      ix = Ax*solp; iy = Ay*solp;
%      Sx = sign(ix); Sy = sign(iy);
      Sx = ix./(sqrt(ix.^2)+epsilon); Sy = iy./(sqrt(iy.^2)+epsilon);
      Ux = epsilon/(ix.^2+epsilon).^(3/2); Uy = epsilon/(iy.^2+epsilon).^(3/2);
      Xr(:,i) = ix; Yr(:,i) = iy;
      
      Sxr(:,i) = Sx; Syr(:,i) = Sy;
      
      % Right hand side
      b = mu*transpose( (Isum.*Sx)'*Ax + (Isum.*Sy)'*Ay );

      Atan = Lhs + mu*( Ax'*(Isum.*Ux)*(Isum.*Ux)'*Ax +...
                        Ay'*(Isum.*Uy)*(Isum.*Ux)'*Ay );
      
      dire = zeros(size(Rhs));
      dire(tokeep) = solp(tokeep) + Atan(tokeep,tokeep)\(b(tokeep)+Rhs(tokeep)); % preconditionned steepest descent direction
%      dire = Lhs*solp + b + Rhs; % steepest descent direction

      sol = solp - relax*dire;
      
      relaxrec(i)  = relax;
      quadratic(i) = .5*sol(tokeep)'*Lhs(tokeep,tokeep)*sol(tokeep) + sol(tokeep)'*Rhs(tokeep);
      absolute(i)  = mu*(Isum'*sqrt((Ax*sol).^2+epsilon)+Isum'*sqrt((Ay*sol).^2+epsilon));
      fctest(i)    = quadratic(i) + absolute(i);
      
      if i>1
         difftest(i) = fctest(i)-fctest(i-1);
      else
         difftest(i) = fctest(i) - fctest0; % First is 0
      end
      
      if i>1 
         crit2(i) = norm(Bxy*(solp-sol))/norm(Bxy*solp);
      else
         crit2(i) = 1;
      end
%      crit(i) = norm((solp-sol))/norm(solp);
      if crit2(i) < 1e-12
         break;
      end
      
%      if i>1
         indexfirst = max(1,i-1);
         criterion = max(difftest(indexfirst:i));
         if criterion <= 0 && relax <= .5
            relax = 2*relax;
         elseif difftest(i) > 0
            relax = .5*relax;
            i = i-1;
            sol = solp; % Cancel the iteration
         end
%      end
      solp = sol;
      i = i+1;
   end
   fctestE = .5*sol'*Lhs*sol + sol'*Rhs + mu*(Isum'*abs(Ax*sol)+Isum'*abs(Ay*sol));
   
elseif normU == 2
   if upper_term == 0
      Mregul = mu*Ax'*Ax + mu*Ay'*Ay;
      sol = zeros(size(Rhs));
      sol(tokeep) = ...
       -(Lhs(tokeep,tokeep) + Mregul(tokeep,tokeep))\Rhs(tokeep);
   else
      sol = -(Lhs + mu*Ax'*Ax + mu*Ay'*Ay)\Rhs;
   end
   
else %normU == 3 : Picard stuff
   L = eye(size(Lhs)); ninfty = 0;
%   L = Bxy'*Bxy; ninfty = 0;
%   L = Ax'*Ax+Ay'*Ay; ninfty = 1;

   [Q,Theta] = eig(Lhs,L);
   Q = Q( 1+ninfty:end , 1+ninfty:end ); Theta = Theta( 1+ninfty:end , 1+ninfty:end );
   Q = Q*(Q'*L( 1+ninfty:end , 1+ninfty:end )*Q)^(-1/2); Q = real(Q);
   thetas = diag(Theta);
   [thetas,Ind] = sort( thetas,'descend' );
   Q = Q(:,Ind);
   Thetas = diag(thetas); 
   Theta = Theta(Ind,Ind);

   figure
   hold on;
   plot(log10(abs(thetas)),'Color','green');
   plot(log10(abs(Q'*Rhs( 1+ninfty:end ))),'Color','red');
   plot(log10(abs((Q'*Rhs( 1+ninfty:end ))./thetas)),'Color','black');
   legend('Singular values','Rhs','sol');

   
   bT = Q'*Rhs( 1+ninfty:end ); bT = bT(1:jmax);
   sol = - Q(:,1:jmax) * (Theta(1:jmax,1:jmax)\bT);

end

%%

% Recover the result
% McCoef  = sol(1:(ordp+1)^2);
% Xx      = sol((ordp+1)^2+1:(ordp+1)^2+(sideL+1)^2);
% Yy      = sol((ordp+1)^2+(sideL+1)^2+1:end);
McCoef = sol;%[ Rhs(1)/Lhs(1,1) ; sol ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if upper_term == 0
   McCoe0 = zeros(size(Rhs));
   McCoe0(tokeep) = -Lhs(tokeep,tokeep)\Rhs(tokeep);
else
   McCoe0 = -Lhs\Rhs; % No regularization
end

% Compare functionnals
no0 = .5*McCoe0'*Lhs*McCoe0 + McCoe0'*Rhs;
re0 = mu*( Isum'*abs(Ax*McCoe0) + Isum'*abs(Ay*McCoe0) );

nor = .5*McCoef'*Lhs*McCoef + McCoef'*Rhs;
reg = mu*( Isum'*abs(Ax*McCoef) + Isum'*abs(Ay*McCoef) );

% Plot the old stuff
nxs = (max(Xs)-min(Xs))/100; nys = (max(Ys)-min(Ys))/100;
X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
solup = zeros(101,101);
solu0 = zeros(101,101);
for k=0:ordp
   for l=0:ordp
      solup = solup + McCoef(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
      solu0 = solu0 + McCoe0(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
   end
end
solup = solup'; % prepare for plot
solu0 = solu0'; %

figure;
hold on;
surf(X,Y,solu0);
shading interp;
colorbar();
axis('equal');

figure;
hold on;
surf(X,Y,solup);
shading interp;
colorbar();
axis('equal');

toc

error0 = norm(solu0(:)-uplo1(:))/norm(uplo1(:));
error1 = norm(solup(:)-uplo1(:))/norm(uplo1(:));
gap    = norm(solup(:)-solu0(:))/norm(solu0(:));

if id_crack == 1 % Identify the crack
   crackID = zeros(size(solup));
   maxto = max(max(solup))*threshold;
   tofind = find( solup>maxto )-1;
   indicei = rem( tofind , size(solup,1) )+1;
   indicej = floor( tofind / size(solup,1) )+1;
   
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

% Plot on the line X = 4
figure;
hold on;
nys = (max(Ys)-min(Ys))/100;
Y = min(Ys):nys:max(Ys); X = pm4;

solup = 0; solu0 = 0;
for k=0:ordp
   for l=0:ordp
      solup = solup + McCoef(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
      solu0 = solu0 + McCoe0(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
   end
end

plot( Y, solup, 'Color', 'black' );
plot( Y, solu0, 'Color', 'blue' );
plot( Y, uplo(3:3:end), 'Color', 'red' );
legend('filtred','unfiltred','reference');
