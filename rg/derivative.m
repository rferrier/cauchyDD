% 03/07/2017 : Correction d'un champ par minimisation de la dérivée
 
clear all;
close all;
 
Norm = 1; ordpD = 2;
 
VAR = load('fields/McCoef1.mat');
McCoef = VAR.McCoef; ordp = VAR.ordp; Xs = VAR.Xs; Ys = VAR.Ys;
X0 = VAR.X0; Y0 = VAR.Y0;
Lx = VAR.Lx; Ly = VAR.Ly;
L1x = VAR.L1x; L2x = VAR.L2x; L1y = VAR.L1y; L2y = VAR.L2y; 

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

%figure;
%hold on;
%surf(X,Y,solup);
%shading interp;
%colorbar();
%axis('equal');

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
            PS0(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) = ...
                 Lx^2*(L2x^ordx - L1x^ordx)/ordx * (L2y^ordy - L1y^ordy)/ordy;
                 % Lx* beacuse there is a variable change x' = Lx*x and y'=Lx*y
         end
      end
   end
end

%    [PassD, ~] = eig(LhsO); % This matrix gives the Legendre basis
%    PassD = eye((ordpD+1)^2);
PassD1 = chol(PS0,'lower');     %CHOLESKI
PassD = inv(PassD1');
 
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
                        Lx^2*k*i*(L2x^(k+i-1)-L1x^(k+i-1))*(L2y^(l+j+1)-L1y^(l+j+1))/...
                                             ((k+i-1)*(l+j+1));
               end
               if l+j>1
                  PS1(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) = ...
                        PS1(j+1+(ordpD+1)*i,l+1+(ordpD+1)*k) +...
                        Ly^2*l*j*(L2x^(k+i+1)-L1x^(k+i+1))*(L2y^(l+j-1)-L1y^(l+j-1))/...
                                             ((k+i+1)*(l+j-1));
               end
            end
         end
      end
   end

   RhsD = zeros((ordpD+1)^2,1); % Will be truncated later : for now, it contains all the coeffs from 0 to ordpD

   LhsD = PassD'*PS1*PassD;
   % We are only using the last vectors from the basis : put zeros
   for k=0:ordpD
      for l=0:ordpD
          if (k<=ordp && l<=ordp)
            LhsD(l+1+(ordpD+1)*k,:) = 0;
            LhsD(:,l+1+(ordpD+1)*k) = 0;
         end
      end
   end

   for k=0:ordpD
      for l=0:ordpD
         for i=0:ordp
            for j=0:ordp

               for ii=0:k
                  for jj=0:l
                     if (k>ordp || l>ordp) && i+ii>1
                        RhsD(l+1+(ordpD+1)*k) = RhsD(l+1+(ordpD+1)*k) - McCoef(1+j+(ordp+1)*i) * ...
                              PassD(jj+1+(ordpD+1)*ii, l+1+(ordpD+1)*k) * ...
                              Lx^2*i*ii*(L2x^(i+ii-1)-L1x^(i+ii-1))*(L2y^(j+jj+1)-L1y^(j+jj+1))/((i+ii-1)*(j+jj+1));
                     end
                     if (k>ordp || l>ordp) && j+jj>1
                        RhsD(l+1+(ordpD+1)*k) = RhsD(l+1+(ordpD+1)*k) - McCoef(1+j+(ordp+1)*i) * ...
                              PassD(jj+1+(ordpD+1)*ii, l+1+(ordpD+1)*k) * ...
                              Ly^2*j*jj*(L2x^(i+ii+1)-L1x^(i+ii+1))*(L2y^(j+jj-1)-L1y^(j+jj-1))/((i+ii+1)*(j+jj-1));
                     end
                  end
               end
               
            end
         end
      end
   end

   % At this point, there are plenty of zero columns (and lines) in LhsD.
   % We are going to add identity on those in order to make the system invertible
   % This has no impact on the solution as the RhsD in front is 0.
   for i=0:ordpD
      for j=0:ordpD
         if norm(LhsD(j+1+(ordpD+1)*i,:)) == 0
            LhsD(j+1+(ordpD+1)*i,j+1+(ordpD+1)*i) = 1;
         end
      end
   end
     
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

%   figure;
%   hold on;
%   surf(X,Y,solupp);
%   shading interp;
%   colorbar();
%   axis('equal');
    
elseif Norm == 1 % L1 norm of the gradient
   % First task : build the matrix giving the derivatives from the coefficients
   nxs = (max(Xs)-min(Xs))/10; nys = (max(Ys)-min(Ys))/10;
   X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
    
   Ax0 = zeros(121,(ordpD+1)^2); Ay0 = zeros(121,(ordpD+1)^2);
   for k=1:ordpD
      for l=0:ordpD
         Amute = Lx*k*(X/Lx-X0)'.^(k-1) * (Y/Lx-Y0).^l;
         Ax0(:,l+1+(ordpD+1)*k) = Amute(:);
         % Lx* beacuse there is a variable change x' = Lx*x and y'=Lx*y
      end
   end
   for k=0:ordpD
      for l=1:ordpD
         Amute = Ly*(X/Lx-X0)'.^k * l*(Y/Lx-Y0).^(l-1);
         Ay0(:,l+1+(ordpD+1)*k) = Amute(:);
         % Lx* beacuse there is a variable change x' = Lx*x and y'=Lx*y
      end
   end
    
   Ax = Ax0*PassD; % Put it in the "Legendre" basis
   Ay = Ay0*PassD;
    
   % ... And re-write over them for the first coefficients (those that are not in the Legendre basis)
   for k=1:ordp
      for l=0:ordp
         Amute = Lx*k*(X/Lx-X0)'.^(k-1) * (Y/Lx-Y0).^l;
         Ax(:,l+1+(ordpD+1)*k) = Amute(:);
      end
   end
   for k=0:ordp
      for l=1:ordp
         Amute = Ly*(X/Lx-X0)'.^k * l*(Y/Lx-Y0).^(l-1);
         Ay(:,l+1+(ordpD+1)*k) = Amute(:);
      end
   end
   
   % Build the matrix of the linear optim problem
   C1 = zeros(4*121,4*(ordpD+1)^2);
   C1( 1:121 , 1:(ordpD+1)^2 )                           = Ax;
   C1( 121+1:2*121 , (ordpD+1)^2+1:2*(ordpD+1)^2 )     = Ax;
   C1( 2*121+1:3*121 , 2*(ordpD+1)^2+1:3*(ordpD+1)^2 ) = Ay;
   C1( 3*121+1:4*121 , 3*(ordpD+1)^2+1:4*(ordpD+1)^2 ) = Ay;
   
   % Do the sum
   toremove = [];
   C = zeros(1,4*(ordpD+1)^2);
   for i=1:4*(ordpD+1)^2
      C(1,i) = sum(C1(:,i));
      if norm(C1(:,i)) == 0
         toremove(end+1) = i;
      end
   end
%   
   % Now, the constraint matrix
   Id = eye((ordpD+1)^2); Im = zeros((ordpD+1)^2);
   Zoo = zeros((ordpD+1)^2); Zno = zeros( 121 , (ordpD+1)^2 );
   Itor = zeros((ordpD+1)^2);
   
%   Itor(toremove,toremove) = eye(size(toremove)); % For the 
   
   for k=0:ordp
      for l=0:ordp
         Im( 1+l+(ordpD+1)*k , 1+l+(ordpD+1)*k ) = 1; % Already-determined coefs
         Itor( 1+l+(ordpD+1)*k , 1+l+(ordpD+1)*k ) = 0; % Remove it from the non-redondant stuff
      end
   end

   Ac = [ Im, -Im, Zoo, Zoo ;...
          Zoo, Zoo, Im, -Im ;...
          Id, -Id, -Id, Id ];
          
   Ac = [Ac;C1];
          
   % Lhs of the constraint
   bc = zeros( 4*121+3*(ordpD+1)^2 , 1 );
   
   a0 = zeros( (ordpD+1)^2 , 1 ); % Already determined coefficients
   for k=0:ordp
      for l=0:ordp
         a0(l+1+(ordpD+1)*k,1) = McCoef(l+1+(ordp+1)*k,1); 
      end
   end
   
   bc( 1:2*(ordpD+1)^2 ) = [a0;a0];
   
   lb = -1e5*ones(4*(ordpD+1)^2,1);
   ub = 1e5*ones(4*(ordpD+1)^2,1); % Bounds
   
   cvect = 'C'*ones( 4*(ordpD+1)^2 , 1 ); cvect = char(cvect);
   ctype1 = 'S'*ones( 3*(ordpD+1)^2 , 1 );
   ctype2 = 'L'*ones( 4*121 , 1 );
   ctype = char( [ctype1;ctype2] );
   
%   ctype = 'U'*ones(4*121 + 4*(ordpD+1)^2,1); ctype = char(ctype);
   % Call the Solver
   param.lpsolver = 1;
   [ sol , minim , errnum , extra ] = glpk ( C' , Ac , bc , lb , ub , ctype , cvect , 1 , param); %(1:27,:)
   
   % Recover the result
   McCoefDt = zeros(4*(ordpD+1)^2,1);
%   McCoefDt(~toremove) = sol;
%   McCoefDt(~toremove) = sol;
   McCoefDt = sol;
   
   McCoefD  = McCoefDt(1:(ordpD+1)^2) - McCoefDt((ordpD+1)^2+1:2*(ordpD+1)^2);
   McCoefD2 = McCoefDt(2*(ordpD+1)^2+1:3*(ordpD+1)^2) - ... 
                     McCoefDt(3*(ordpD+1)^2+1:4*(ordpD+1)^2);
                     
   % Actually, we should have McCoefD = McCoefD2
   
   solupp = solup;
   
else % No gradient stuff
   solupp = solup;
   McCoefD = zeros((ordpD+1)^2,1);
end
 
toc
 
 % Plot on the line X = 4
%figure;
%hold on;
%nys = (max(Ys)-min(Ys))/100;
%Y = min(Ys):nys:max(Ys); X = 4;
%
%solup = 0;
%for k=0:ordp
%   for l=0:ordp
%      solup = solup + McCoef(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
%   end
%end
%
%plot( Y, solup, 'Color', 'black' );
%
%solupp = solup;
%for k=0:ordpD
%   for l=0:ordpD
%
%      for ii=0:k
%         for jj=0:l
%            solupp = solupp + PassD(jj+1+(ordpD+1)*ii, l+1+(ordpD+1)*k) * ...
%                    McCoefD(1+l+(ordpD+1)*k) .* (X/Lx-X0)'.^ii * (Y/Lx-Y0).^jj;
%         end
%      end
%  end
%end
%plot( Y, solupp, 'Color', 'green' );
