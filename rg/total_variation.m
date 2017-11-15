% 13/11/2017 : Identification par minimisation d'une fonctionnelle
% régularisée par Total Variation
 
clear all;
close all;

sideL = 20; % Nb of points on the domain (for integral approximation)
pm4   = -4; % Position of the 1d line
mu    = 1e-2;%5e-4;%1e-4;%5e-5; % Regularization parameter
ordp  = 8;

VAR = load('fields/AnBm15b10.mat');
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

for k=1:ordp
   for l=0:ordp
      Amute = 1/Lx*k*(X/Lx-X0)'.^(k-1) * (Y/Lx-Y0).^l;
      Ax(:,l+1+(ordp+1)*k) = Amute(:);
      Bmute = (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
      Bxy(:,l+1+(ordp+1)*k) = Bmute(:);
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

% Build the matrix of the linear optim problem
% M = [ zeros(1,(ordp+1)^2) , Isum' , Isum' ]'; %ones(1,2*(sideL+1)^2)

% Now, the constraint matrix
Id = eye((sideL+1)^2); Im = zeros(0,(ordp+1)^2);
Zoo = zeros((sideL+1)^2); Zno = zeros( (sideL+1)^2 , (ordp+1)^2 );
   
%   Itor(toremove,toremove) = eye(size(toremove)); % For the 
ind = 1;
for k=0:ordp
   for l=0:ordp
      Im( ind , 1+l+(ordp+1)*k ) = 1; % Already-determined coefs
      ind = ind+1;
   end
end

Zoom = zeros(size(Im,1),(sideL+1)^2);
   
C =  [ Ax, -Id, Zoo ;...
       -Ax, -Id, Zoo ;...
       Ay, Zoo, -Id ;...
       -Ay, Zoo, -Id];
          
% Lhs of the constraint
d = zeros( 4*(sideL+1)^2 , 1 );
  
% Regular minimization matrix
Ire = norm(Lhs)*0 * eye(size(Zoo)); % To regularize the linear inversions
A  = [ Lhs, Zno', Zno' ;...
       Zno, Ire, Zoo ;...
       Zno, Zoo, Ire ];
% Am1 = pinv(A);
b   = [ -Rhs ; mu*Isum ; mu*Isum ];
bt1 = mu*[ zeros(size(Rhs)) ; Isum ; Isum ];
bt2 = [ -Rhs ; zeros(size(Isum)) ; zeros(size(Isum)) ];

ndof = size(A,1); nimp = size(C,1);

%% Uzawa trial
% sol  = Am1*b;
% % % kernel postpro
% % Xx = abs(Ax*sol(1:(ordp+1)^2));
% % Yy = abs(Ay*sol(1:(ordp+1)^2));
% % sol((ordp+1)^2+1:(ordp+1)^2+(sideL+1)^2) = Xx;
% % sol((ordp+1)^2+(sideL+1)^2+1:end)        = Yy;
% 
% f = zeros(nimp,1);
% 
% % Use Uzawa to solve the system
% k = 1e-3;% norm(C,'fro');
% crit = zeros(100,1); critf = zeros(100,1);
% for i=1:100
%    delta = C*sol-d;
% 
%    fp = f;
%    for j=1:nimp
%       if f(j)-k*delta(j) < 0
%          f(j) = f(j)-k*delta(j); % Negative part
%       else
%          f(j) = 0;
%       end
%    end
% 
%    solp = sol(1:(ordp+1)^2); solp0 = sol;
%    sol  = Am1 * (b + C'*f);
%    
% %    % kernel postpro
% %    Xx = abs(Ax*sol(1:(ordp+1)^2));
% %    Yy = abs(Ay*sol(1:(ordp+1)^2));
% %    sol((ordp+1)^2+1:(ordp+1)^2+(sideL+1)^2) = Xx;
% %    sol((ordp+1)^2+(sideL+1)^2+1:end)        = Yy;
% 
%    critf(i) = norm(fp-f)/norm(fp);
%    crit(i)  = norm(solp-sol(1:(ordp+1)^2))/norm(solp);
%    bug
%    if crit(i) < 1e-12
%      break;
%   end
% end

%% Homemade stuff

sol   = zeros(size(Rhs)); % Initialize
crit  = zeros(100,1);
relax = .25;
relaxrec = zeros(100,1); % Store the relax parameter (debug)

for i=1:100
   % Compute signs
   Sx = sign(Ax*sol); Sy = sign(Ay*sol);
   
   % Right hand side
   b = mu*transpose( (Isum.*Sx)'*Ax + (Isum.*Sy)'*Ay );
   
   % Solve again
   solp = sol;
   sol = Lhs\(-b-Rhs);
   
   if i>1 sol = relax*sol + (1-relax)*solp; end % Relaxation step
   relaxrec(i) = relax;
   
   crit(i) = norm(Bxy*(solp-sol))/norm(Bxy*solp);
   if crit(i) < 1e-12
      break;
   end
   
   if i>2
      olddiff = crit(i-1)-crit(i-2);
      newdiff = crit(i)-crit(i-1);
      if crit(i) < crit(i-1) && abs(olddiff) < abs(newdiff)
         relax = 2*relax;
      else
         relax = .5*relax;
      end
   end
%    relax
end

%%

% Recover the result
% McCoef  = sol(1:(ordp+1)^2);
% Xx      = sol((ordp+1)^2+1:(ordp+1)^2+(sideL+1)^2);
% Yy      = sol((ordp+1)^2+(sideL+1)^2+1:end);
McCoef = sol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
McCoe0 = -Lhs\Rhs; % No regularization

% Compare functionnals
no0 = McCoe0'*Lhs*McCoe0 + McCoe0'*Rhs;
re0 = mu*( Isum'*abs(Ax*McCoe0) + Isum'*abs(Ay*McCoe0) );

nor = McCoef'*Lhs*McCoef + McCoef'*Rhs;
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

error0 = norm(solu0-uplo1)/norm(uplo1);
error1 = norm(solup-uplo1)/norm(uplo1);
gap    = norm(solup-solu0)/norm(solu0);

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
