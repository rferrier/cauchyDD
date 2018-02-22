% 30/11/2017 : Identification par minimisation d'une fonctionnelle
% avec valeurs imposées sur les bords, et positivité
 
clear all;
close all;

sideL      = 10; % Nb of points on the domain (for integral approximation)
pm4        = 4; % Position of the 1d line
kpen       = 0;%Penalization parameter
ordp       = 13;
upper_term = 0;

%VAR = load('fields/AnBnoplane615.mat');
%VAR = load('fields/AnBm15b2.mat');
VAR = load('fields/AnBm15.mat');
%VAR = load('fields/AnBs15.mat');
%VAR = load('fields/AnBs15b1.mat');
solpoly = csvread('fields/rg3d_four.csv');
linesolpoly = load('fields/fourier3D2D.mat');
lpoly = linesolpoly.solu;

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
caxis( [-0.004497, 0.020062] );
colorbar('SouthOutside');

tic

% First task : build the matrix giving the values from the coefficients
nxs = (max(Xs)-min(Xs))/sideL; nys = (max(Ys)-min(Ys))/sideL;
X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
mX = min(Xs)*ones(1,sideL+1); MX = max(Xs)*ones(1,sideL+1);
mY = min(Ys)*ones(1,sideL+1); MY = max(Ys)*ones(1,sideL+1);

Xb = [X,X,mX,MX]; Yb = [mY,MY,Y,Y];
Axy = zeros((sideL+1)^2,(ordp+1)^2); % Value
Bxy = zeros(4*(sideL+1),(ordp+1)^2); % Value on the bound

for k=0:ordp
   for l=0:ordp
      Amute = (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
      Bmute = (Xb/Lx-X0).^k .* (Yb/Lx-Y0).^l;
      Axy(:,l+1+(ordp+1)*k) = Amute(:);
      Bxy(:,l+1+(ordp+1)*k) = Bmute;
   end
end

% Reference stuff
McCoe0 = zeros(size(Rhs));
McCoe0(tokeep) = -Lhs(tokeep,tokeep)\Rhs(tokeep);

% Boundary stuff
Cbound = Bxy'*Bxy;
McCoef = zeros(size(Rhs));
McCoef(tokeep) = -( Lhs(tokeep,tokeep) + kpen*Cbound(tokeep,tokeep) )\Rhs(tokeep);

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
caxis( [-0.004497, 0.020062] );
colorbar('SouthOutside');

figure;
hold on;
surf(X,Y,solup);
shading interp;
colorbar();
axis('equal');
caxis( [-0.004497, 0.020062] );
colorbar('SouthOutside');

figure;
hold on;
surf(X,Y,solpoly);
shading interp;
colorbar();
axis('equal');
caxis( [-0.004497, 0.020062] );
colorbar('SouthOutside');

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
plot( Y, lpoly, 'Color', 'green' );
plot( Y, uplo(3:3:end), 'Color', 'red' );
legend('filtred','unfiltred','fourier','reference');
