% 03/07/2017 : Correction d'un champ par minimisation de la dérivée
 
clear all;
close all;
 
Norm = 2;
ordpD = 12;
sideL = 10; % Nb of points on the domain (for integral approximation)
pm4   = 4; % Position of the 1d line
 
%VAR = load('fields/McCoef5.mat');
%VAR = load('fields/McCoefs7.mat');
%VAR = load('fields/McCoef8r.mat');
%VAR = load('fields/McCoef7r6.mat');
%VAR = load('fields/McCoef7rb1.mat');
VAR = load('fields/McCoef6rb1.mat');
McCoef = VAR.McCoef; ordp = VAR.ordp; Xs = VAR.Xs; Ys = VAR.Ys;
X0 = VAR.X0; Y0 = VAR.Y0;
Lx = VAR.Lx; Ly = VAR.Ly;
L1x = VAR.L1x; L2x = VAR.L2x; L1y = VAR.L1y; L2y = VAR.L2y;

%REF = load('fields/ref_mult_u.mat');
REF = load('fields/ref_mult.mat');
%REF = load('fields/ref_mults.mat');
%REF = load('fields/ref_mult6.mat');
uplo1 = REF.uplo1; uplo = REF.uplo;

% Plot the old stuff
nxs = (max(Xs)-min(Xs))/100; nys = (max(Ys)-min(Ys))/100;
X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
solup = zeros(101,101);
for k=0:ordp
   for l=0:ordp
      solup = solup + McCoef(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
   end
end
solup = solup'; % prepare for plot

figure;
hold on;
surf(X,Y,solup);
shading interp;
colorbar();
axis('equal');

% The reference
figure;
hold on;
surf(X,Y,uplo1);
shading interp;
colorbar();
axis('equal');

error0 = norm( solup-uplo1 ) / norm(uplo1);

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

   LhsD = PassD'*PS1*PassD;
   RhsD1 = PassD'*RhsR;
   % We are only using the last vectors from the basis : put zeros
   for k=0:ordpD
      for l=0:ordpD
          if (k<=ordp && l<=ordp)
            LhsD(l+1+(ordpD+1)*k,:) = 0;
            LhsD(:,l+1+(ordpD+1)*k) = 0;
         end
      end
   end

%   for k=0:ordpD
%      for l=0:ordpD
%         for i=0:ordp
%            for j=0:ordp
%
%               for ii=0:k
%                  for jj=0:l
%                     if (k>ordp || l>ordp) && i+ii>1
%                        RhsD(l+1+(ordpD+1)*k) = RhsD(l+1+(ordpD+1)*k) - McCoef(1+j+(ordp+1)*i) * ...
%                              PassD(jj+1+(ordpD+1)*ii, l+1+(ordpD+1)*k) * ...
%                              1/Lx^2*i*ii*(L2x^(i+ii-1)-L1x^(i+ii-1))*(L2y^(j+jj+1)-L1y^(j+jj+1))/((i+ii-1)*(j+jj+1));
%                     end
%                     if (k>ordp || l>ordp) && j+jj>1
%                        RhsD(l+1+(ordpD+1)*k) = RhsD(l+1+(ordpD+1)*k) - McCoef(1+j+(ordp+1)*i) * ...
%                              PassD(jj+1+(ordpD+1)*ii, l+1+(ordpD+1)*k) * ...
%                              1/Ly^2*j*jj*(L2x^(i+ii+1)-L1x^(i+ii+1))*(L2y^(j+jj-1)-L1y^(j+jj-1))/((i+ii+1)*(j+jj-1));
%                     end
%                  end
%               end
%               
%            end
%         end
%      end
%   end

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

   % Add it to solup
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
   
   error = norm( solupp-uplo1 ) / norm(uplo1);

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
   Id = eye((sideL+1)^2); Im = zeros(0,(ordpD+1)^2);
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
   
   % Call the Solver
   param.lpsolver = 1;
   [ sol , minim , errnum , extra ] = glpk ( C' , Ac , bc , lb , ub , ctype , cvect , 1 , param); %(1:27,:)
   
   % Recover the result
%   McCoefDt = zeros(4*(ordpD+1)^2,1);
%   McCoefDt(~toremove) = sol;
%   McCoefDt(~toremove) = sol;

   McCoefDt  = sol(1:(ordpD+1)^2);
   Xx        = sol((ordpD+1)^2+1:(ordpD+1)^2+(sideL+1)^2);
   Yy        = sol((ordpD+1)^2+(sideL+1)^2+1:end);
   
   % Remove the already existing coeffs (as it will be added later)
   McCoefD = McCoefDt;
   McCoef0 = McCoefDt-McCoefDt;
   for k=0:ordp
      for l=0:ordp
         McCoefD(l+1+(ordpD+1)*k,1) = 0;
         McCoef0(l+1+(ordpD+1)*k,1) = McCoef(l+1+(ordp+1)*k,1);
      end
   end
   
   % Debug : see the value without any McCoefDt
   Xx0 = abs(Ax*McCoef0); Yy0 = abs(Ay*McCoef0);
   val0 = C*[McCoef0;Xx0;Yy0];
   
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
   
   error = norm( solupp-uplo1 ) / norm(uplo1);
   
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

plot( Y, solup, 'Color', 'black' );

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
