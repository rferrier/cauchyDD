% 13/11/2017 : Identification par minimisation d'une fonctionnelle
% régularisée par Total Variation Latin dual
 
clear all;
close all;

sideL      = 10; % Nb of points on the domain (for integral approximation)
pm4        = 4; % Position of the 1d line
%mu         = 1e-3;%1e-4;%1e-3;%.7e-2;%5e-4;%1e-4;%5e-5; % Regularization parameter
mu         = 1e-1;
ordp       = 6;
normU      = 2;   % Use L1 or L2 minimization
upper_term = 0;
klatin     = 0;%1e2; % Latin Research direction (0 means norm(GAG,'fro')/10)
nlatin     = 500;
omega      = .5; % Relaxation parameter
id_crack   = 0;  % Should I identify the crack (binary stuff) ?
threshold  = .1;  % Threshold for the crack

%VAR = load('fields/AnBm15b2.mat');
%VAR = load('fields/AnBm15.mat');
%VAR = load('fields/AnBs15.mat');
%VAR = load('fields/AnBs15b1.mat');
VAR = load('fields/AnBnoplane2615.mat');
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

beginning = min(Xs); ending = max(Xs);
%beginning = -1.0375; ending = 6.5547;
nxs = (ending-beginning)/100; nys = (max(Ys)-min(Ys))/100;
X = beginning:nxs:ending; Y = min(Ys):nys:max(Ys);
% The reference
figure;
hold on;
surf(X,Y,uplo1);
shading interp;
colorbar();
axis('equal');
%caxis( [-0.011621, 0.016879] );
%caxis( [-8.9798e-04, 0.0035356 ]);
colorbar('SouthOutside');

% Orthogonality is useless in this case
tic

% First task : build the matrix giving the derivatives from the coefficients
nxs = (ending-beginning)/sideL; nys = (max(Ys)-min(Ys))/sideL;
X = beginning:nxs:ending; Y = min(Ys):nys:max(Ys);
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

% Build the matrix of the linear optim problem
% M = [ zeros(1,(ordp+1)^2) , Isum' , Isum' ]'; %ones(1,2*(sideL+1)^2)

% Latin for the dual problem
if normU == 1 && mu ~= 0
   M = [ zeros(1,(ordp+1)^2) , Isum' , Isum' ]';
   
   if upper_term==0
      G = [Ax(:,tokeep);Ay(:,tokeep)];
      invA = Lhs(tokeep,tokeep)\eye(size(tokeep,2)); % This is better than inv for some reason
      b = -Rhs(tokeep);
   else
      G = [Ax;Ay];
      invA = Lhs\eye(size(Lhs)); % This is better than inv for some reason
      b = -Rhs;
   end
   
   GAG   = G*invA*G';
   GAbmu = G*invA*b/mu;

   sol0 = zeros(size(GAG,1),1);
   %options.MaxIter = 1000;
   
%   [ sol2 , obj , info , lambda ] = qp (sol0, GAG, GAbmu, Akk, bkk, lb, ub, options);
   
   % Or Image Uzawa method
   K = null(GAG); P = K*((K'*K)\K'); %(normally, (K'*K) = 1 but never trust the matrix)
   sol2 = GAG\GAbmu; sol2 = sol2-P*sol2; % Yes, there is a warn here, but don't care, I manage the stuff
   fsol = zeros(size(sol2));

   if klatin == 0 % Automatic choose
      k = norm(GAG,'fro')/10;
   else
      k = klatin;
   end

   C = [ eye(size(GAG,1)) ; -eye(size(GAG,1)) ];
   nimp = size(C,1);    % Nb of imposed relations
   ndof = size(sol2,1); % Nb of dofs
   
   IsumIsum = [Isum;Isum];% - P*[Isum;Isum];
   d = [IsumIsum;IsumIsum];
   f = zeros(2*size(sol0,1),1); fc = zeros(size(C,1),1);
   critd = [0;0];

   stag = zeros(nlatin,1);

   % Solve min .5*w'*GAG*w - w'*GAbmu, with |w| <= IsumIsum with the LaTIn method
   for i=1:nlatin
      %% Local step (dof per dof status algo)
      % First trial : fhat needs to be < 0
      uhat = IsumIsum;
      fhat = fsol + k*(uhat - sol2);
      tochange = find(fhat > 0);
      % Second trial : fhat needs to be > 0
      uhat(tochange) = -IsumIsum(tochange);
      fhat(tochange) = fsol(tochange) + k*(uhat(tochange) - sol2(tochange));
      tochange2 = find(fhat(tochange) < 0);
      % Last solution : uhat is inbetween
      fhat(tochange(tochange2)) = 0;
      uhat(tochange(tochange2)) = sol2(tochange(tochange2)) - k*fsol(tochange(tochange2));

      %% Linear step
      sol2t = (GAG+k*eye(size(GAG))) \ (GAbmu + fhat + k*uhat);
      
      % Evaluate stagnation criterion
      stag(i) = norm(sol2-sol2t)/norm(sol2);

      % Relaxation step
      sol2 = omega*sol2t + (1-omega)*sol2;
   end
   
   sol3  = sol2;

   if upper_term == 0
      sol = zeros(size(Lhs,1),1);
      sol(tokeep) = Lhs(tokeep,tokeep)\(b-mu*G'*sol3);
   else
      sol = Lhs\(b-mu*G'*sol3);
   end
   
   fctestE = .5*sol(tokeep)'*Lhs(tokeep,tokeep)*sol(tokeep) +...
             sol(tokeep)'*Rhs(tokeep) +...
             mu*(Isum'*abs(Ax(:,tokeep)*sol(tokeep))+...
                 Isum'*abs(Ay(:,tokeep)*sol(tokeep)));
                 
   nor = .5*sol'*Lhs*sol+ sol'*Rhs;
   reg = mu*(Isum'*abs(Ax*sol)+Isum'*abs(Ay*sol));


else %if normU == 2
   if upper_term == 0
      Mregul = mu*Ax'*Ax + mu*Ay'*Ay;
      sol = zeros(size(Rhs));
      sol(tokeep) = ...
       -(Lhs(tokeep,tokeep) + Mregul(tokeep,tokeep))\Rhs(tokeep);
   else
      sol = -(Lhs + mu*Ax'*Ax + mu*Ay'*Ay)\Rhs;
   end
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
nxs = (ending-beginning)/100; nys = (max(Ys)-min(Ys))/100;
X = beginning:nxs:ending; Y = min(Ys):nys:max(Ys);
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
surf(X,Y,solu0);%surf(X,Y,solu0(:,end:-1:1));
shading interp;
colorbar();
axis('equal');
%caxis( [-0.011621, 0.016879] );
caxis( [-0.0031520, 0.0044433 ]);
colorbar('SouthOutside');

figure;
hold on;
surf(X,Y,solup);%(:,end:-1:1)); % /!\ May be necessary to revert that axis
shading interp;
colorbar();
axis('equal');
%caxis( [-0.011621, 0.016879] );
caxis( [-0.0031520, 0.0044433 ]);
colorbar('SouthOutside');

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
