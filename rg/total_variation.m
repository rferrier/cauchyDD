% 13/11/2017 : Identification par minimisation d'une fonctionnelle
% régularisée par Total Variation
 
clear all;
close all;

sideL      = 10; % Nb of points on the domain (for integral approximation)
pm4        = -4; % Position of the 1d line
mu         = 1e-3;%1e-4;%1e-3;%.7e-2;%5e-4;%1e-4;%5e-5; % Regularization parameter
%mu         = 1e-1;
ordp       = 12;
normU      = 1;   % Use L1 or L2 minimization
jmax       = 20;  % Coefficient for Picard stuff (if normU == 3)
upper_term = 0;
kUzawa     = 1e-1;%2e-7; % Uzawa parameter
nuzawa     = 5000;
id_crack   = 0;  % Should I identify the crack (binary stuff) ?
threshold  = .1;  % Threshold for the crack

VAR = load('fields/AnBm15b2.mat');
%VAR = load('fields/AnBm15.mat');
%VAR = load('fields/AnBs15.mat');
%VAR = load('fields/AnBs15b1.mat');
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
Is2 = ones(sideL+1,1); %Is2(1) = .5; Is2(end) = .5;
Isum = Is2*Is2'; Isum = Isum(:);
%Isum = Isum*nxs*nys;

% Build the matrix of the linear optim problem
% M = [ zeros(1,(ordp+1)^2) , Isum' , Isum' ]'; %ones(1,2*(sideL+1)^2)

%% Now, the constraint matrix
%Id = eye((sideL+1)^2); Im = zeros(0,(ordp+1)^2);
%Zoo = zeros((sideL+1)^2); Zno = zeros( (sideL+1)^2 , (ordp+1)^2 );
%   
%%   Itor(toremove,toremove) = eye(size(toremove)); % For the 
%ind = 1;
%for k=0:ordp
%   for l=0:ordp
%      Im( ind , 1+l+(ordp+1)*k ) = 1; % Already-determined coefs
%      ind = ind+1;
%   end
%end

%Zoom = zeros(size(Im,1),(sideL+1)^2);
   
%C =  [ Ax, -Id, Zoo ;...
%       -Ax, -Id, Zoo ;...
%       Ay, Zoo, -Id ;...
%       -Ay, Zoo, -Id];
          
% Lhs of the constraint
%d = zeros( 4*(sideL+1)^2 , 1 );
  
% Regular minimization matrix
%Ire = Zoo;%norm(Lhs)*0 * eye(size(Zoo)); % To regularize the linear inversions
%A  = [ Lhs, Zno', Zno' ;...
%       Zno, Ire, Zoo ;...
%       Zno, Zoo, Ire ];

% Am1 = pinv(A);
%b   = [ -Rhs ; mu*Isum ; mu*Isum ];
%bt1 = mu*[ zeros(size(Rhs)) ; Isum ; Isum ];
%bt2 = [ -Rhs ; zeros(size(Isum)) ; zeros(size(Isum)) ];

%ndof = size(A,1); nimp = size(C,1);

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

% Quadratic programming / Uzawa for the dual problem
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
   sol2 = GAG\GAbmu; sol2 = sol2-P*sol2;
   k = kUzawa;
   C = [ eye(size(GAG,1)) ; -eye(size(GAG,1)) ]; 
   
   IsumIsum = [Isum;Isum];% - P*[Isum;Isum];
   d = [IsumIsum;IsumIsum];
   f = zeros(2*size(sol0,1),1); fc = zeros(size(C,1),1);
   critd = [0;0];
   for i=1:nuzawa
      delta = C*sol2-d;

      fp = f;
      fn = fc-k*delta;

      fc = .5*( fn - abs(fn) ); f = C'*fc; %Negative part
      critd(i) = max(delta+abs(delta));

      solp = sol2;
      sol2 = GAG\(GAbmu+f); % pinv is used in this case
      sol2 = sol2-P*sol2;
   end
   
   if upper_term == 0
      sol = zeros(size(Lhs,1),1);
      sol(tokeep) = Lhs(tokeep,tokeep)\(b-mu*G'*sol2);
   else
      sol = Lhs\(b-mu*G'*sol2);
   end
   
   fctestE = .5*sol(tokeep)'*Lhs(tokeep,tokeep)*sol(tokeep) +...
             sol(tokeep)'*Rhs(tokeep) +...
             mu*(Isum'*abs(Ax(:,tokeep)*sol(tokeep))+...
                 Isum'*abs(Ay(:,tokeep)*sol(tokeep)));
                 
   nor = .5*sol'*Lhs*sol+ sol'*Rhs;
   reg = mu*(Isum'*abs(Ax*sol)+Isum'*abs(Ay*sol));
% Quadratic programming
elseif normU == 11
    M = [ zeros(1,(ordp+1)^2) , Isum' , Isum' ]'; %ones(1,2*(sideL+1)^2)
   
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
   d   = zeros( 4*(sideL+1)^2 , 1 );
   dkk = -1e5*ones(4*(sideL+1)^2,1);
     
   % Regular minimization matrix
%   Ire = Zoo;%norm(Lhs)*0 * eye(size(Zoo)); % To regularize the linear inversions
   A  = [ Lhs, Zno', Zno' ;...
          Zno, Zoo, Zoo ;...
          Zno, Zoo, Zoo ];
   b   = [ -Rhs ; mu*Isum ; mu*Isum ];
   
   % Hack zero therms
   if upper_term == 0
      Akk = zeros(  size(toremove,2), 2*(sideL+1)^2+(ordp+1)^2 );
      
      for i=1:size(toremove,2)
         Akk(i,toremove(i)) = 1;
      end
      
      bkk = zeros(size(toremove,2),1);
   else
      Akk = zeros(0,2*(sideL+1)^2+(ordp+1)^2);
      bkk = zeros(0,1);
   end
     
   lb = -1e5*ones(2*(sideL+1)^2+(ordp+1)^2,1);
   ub = 1e5*ones(2*(sideL+1)^2+(ordp+1)^2,1); % Bounds

   if upper_term == 0
      McCoe0 = zeros(size(Rhs));
      McCoe0(tokeep) = -Lhs(tokeep,tokeep)\Rhs(tokeep);
   else
      McCoe0 = -Lhs\Rhs; % No regularization
   end
   Xx0 = abs(Ax*McCoe0); Yy0 = abs(Ay*McCoe0);
   sol0 = [ McCoe0 ; Xx0 ; Yy0 ];
   
   [ sol2 , obj , info , lambda ] = qp (sol0, A, b, Akk, bkk, lb, ub, dkk, C, d);
   sol = sol2(1:(ordp+1)^2);

% Uzawa with Chambolle dual whatsoever stuff
elseif normU == 1 && mu > 0
   
   % First : the constraint
   G = Bxy;
   C = [G;-G]; nimp = size(C,1); ndof = size(C,2);
   d = [Isum;Isum];%ones(nimp,1);
   A = Lhs;
   b = -Rhs;

   sol = zeros(ndof,1); sol0 = sol; solD = sol;
   
   solD(tokeep) = b(tokeep)/mu; sol0 = solD;
   f   = zeros(ndof,1); fc = zeros(nimp,1);
   
   lam = mu;

   % C*sol - d <= 0;
%   Cp = pinv(C); %k = Cp'*A*Cp/2;% Optimal stiffness
   k = 5e-3;
   
   crit = zeros(2,1); critf = zeros(2,1);
   for i=1:1000000
      delta = C*solD-d;
    
      fp = f;
%      fn = f-k*sign(delta).*sqrt(abs(delta));
      fn = fc-k*delta;

      fc = .5*( fn - abs(fn) ); f = C'*fc; %Negative part
    
      solp          = solD;
      solD(tokeep)  = (b(tokeep)/lam + f(tokeep));
    
%      if norm(fp)>0
         critf(i) = norm(fp-f)/norm(fp);
%      end
      crit(i)  = norm(solp-solD)/norm(solp);
      critd(i) = max(delta+abs(delta)); % Quite time-consuming
   %   bug
      if critd(i) == 0
        break;
      end
      if critd(i)/critd(1) < 1e-9
        break;
      end
      
      % K-adaptativity
%      if i>3
%        if critd(i)<critd(i-1) && critd(i-1)<critd(i-2) && critd(i-2)<critd(i-3)
%           k = 2*k;
%        elseif critd(i)>critd(i-1)
%           k = k/10;
%        end
%     end

%      if i>1
%        if critd(i) > critd(i-1)
%           k = k/2;
%        elseif critd(i) > 3/4*critd(i-1) && critd(i) < 7/8*critd(i-1)
%           k = 1.5*k;
%        end
%     end
%     k
   end
   critd = critd/critd(1);
%   bug;
   % Re-dualize
%   solD(1) = 0; % Because of the Kernel of G
   sol = zeros(ndof,1);
   sol(tokeep) = A(tokeep,tokeep)\(b(tokeep) - lam*solD(tokeep));
%   sol = lam*solD;
   fctestE = .5*sol'*Lhs'*Lhs*sol + sol'*Lhs'*Rhs +...
                  mu*(Isum'*abs(Ax*sol)+Isum'*abs(Ay*sol));
   
   sol0 = A\b;   
   fctest0 = .5*sol0'*Lhs'*Lhs*sol0 + sol0'*Lhs'*Rhs +...
                  mu*(Isum'*abs(Ax*sol0)+Isum'*abs(Ay*sol0));

%elseif normU == 11 && mu > 0 % Fucking bugged stuff
%   
%   % First : the constraint
%   G = [Ax;Ay];
%   C = [G;-G]; nimp = size(C,1); ndof = size(C,2);
%   d = ones(nimp,1);
%   A = Lhs;
%   b = -Rhs;
%   
%%   gg = G\G;
%   
%%   imGpG = eye(ndof); imGpG(1,1) = 0; % Build the Image of G (constant term)
%   td = [1:ndof];
%   sol = zeros(ndof,1); sol0 = sol; solD = sol;
%   
%%   At = imGpG*A*imGpG;
%   At = A(td,td);
%   
%%   Am1 = inv(A); % This one is small
%   Amp1 = At\eye(size(A(td,td))); % Safe way to invert
%%   bug
%   solD(td) = Amp1 * (b(td)/mu); sol0 = solD;
%   f   = zeros(ndof,1); fc = zeros(nimp,1);
%   
%   lam = mu;
%
%   % C*sol - d <= 0;
%   Cp = pinv(C); %k = Cp'*A*Cp/2;% Optimal stiffness
%   k = 3e-5;
%   
%   crit = zeros(2,1); critf = zeros(2,1);
%   for i=1:1000
%      delta = C*solD-d;
%    
%      fp = f;
%
%%      fn = f-k*sign(delta).*sqrt(abs(delta));
%      fn = fc-k*delta;
%%      fn = fc- C'\(A*(C\delta))/10;  % Optimal stifness
%      fc = .5*( fn - abs(fn) ); f = C'*fc;
%    
%      solp = solD;
%      solD(td)  = Amp1 * (b(td)/lam + f(td));
%    
%      critf(i) = norm(fp-f)/norm(fp);
%      crit(i)  = norm(solp-solD)/norm(solp);
%      critd(i) = max(delta+abs(delta)); % Quite time-consuming
%   %   bug
%      if crit(i) < 1e-12
%        break;
%      end
%   end
%%   bug;
%   % Re-dualize
%%   solD(1) = 0; % Because of the Kernel of G
%   sol = A\b - lam*solD;
%%   sol = lam*solD;
%   fctestE = .5*sol'*Lhs*sol + sol'*Rhs +...
%                  mu*(Isum'*abs(Ax*sol)+Isum'*abs(Ay*sol));

%   nxs = (max(Xs)-min(Xs))/100; nys = (max(Ys)-min(Ys))/100;
%   X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
%   soluD = zeros(101,101);
%   for k=0:ordp
%      for l=0:ordp
%         soluD = soluD + solD(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
%      end
%   end
%   soluD = soluD'; % prepare for plot
%   figure;
%   hold on;
%   surf(X,Y,soluD);
%   shading interp;
%   colorbar();
%   axis('equal');
                  
%% Homemade stuff (deactivated)
elseif normU == 1
   sol   = zeros(size(Rhs),1); % Initialize
   crit  = zeros(100,1); crit2 = zeros(100,1); diff = zeros(100,1);
   relax = 1;
   sol0 = -Lhs\Rhs;
   
%   lim  = 1e-2; % Regularization
%   limx = lim*Isum'*abs(Ax*sol0)/sideL^2/nxs/nys; % Limit for each term in Ax
%   limy = lim*Isum'*abs(Ay*sol0)/sideL^2/nxs/nys; % Limit for each term in Ay
   
   Sxr = zeros( size(Ax*sol), 100); Syr = zeros( size(Ay*sol), 100);
   Xr  = zeros( size(Ax*sol), 100); Yr  = zeros( size(Ay*sol), 100); 
   relaxrec = zeros(100,1); % Store the relax parameter (debug)
   fctest = zeros(100,1); difftest = zeros(100,1);
   
   sol = sol0; 
   solp = sol;
   
%   fctest0 = .5*sol'*Lhs'*Lhs*sol + sol'*Lhs'*Rhs + .5*Rhs'*Rhs +...
%                  mu*(Isum'*abs(Ax*sol)+Isum'*abs(Ay*sol));
                  
   fctest0 = .5*sol'*Lhs*sol + sol'*Rhs +...
                  mu*(Isum'*abs(Ax*sol)+Isum'*abs(Ay*sol));
   
   i = 1;
   for j=1:1000
      % Compute signs
%      Sxr(:,i) = Sx; Syr(:,i) = Sy;
      ix = Ax*solp; iy = Ay*solp;
      Sx = sign(ix); Sy = sign(iy);
      Xr(:,i) = ix; Yr(:,i) = iy;
      
%      % Small values correction (regularization)
%%      mSx = mean(abs(Ax*solp)); mSy = mean(abs(Ay*sol));
%      toZeroIfx = find(abs(ix) < limx);%lim*mSx);
%      toZeroIfy = find(abs(iy) < limy);%lim*mSy);
%      Sx(toZeroIfx) = ix(toZeroIfx)/limx; Sy(toZeroIfy) = iy(toZeroIfy)/limy;
      
      Sxr(:,i) = Sx; Syr(:,i) = Sy;
      
      % Right hand side
      b = mu*transpose( (Isum.*Sx)'*Ax + (Isum.*Sy)'*Ay );
      
      % Solve again
%      sol  = Lhs\(-b-Rhs);

      dire = solp + Lhs\(b+Rhs); % preconditionned steepest descent direction
%      dire = Lhs*solp + b + Rhs; % steepest descent direction
      
%      % Randomize
%      bcand = zeros( size(sol,1) , 10*size(sol,1)+1 ); socand = zeros(size(bcand));
%      fcand = zeros(1,10*size(sol,1)+1);
%      for k=1:10*size(sol,1)
%         bcand(:,k)  = abs( rand(size(sol,1),1) ) .* dire;
%         socand(:,k) = solp - relax*bcand(:,k);
%         fcand(k) = .5*socand(:,k)'*Lhs*socand(:,k) + socand(:,k)'*Rhs + ...
%                    mu*(Isum'*abs(Ax*socand(:,k))+Isum'*abs(Ay*socand(:,k)));
%      end
%      bcand(:,end) = dire;
%      socand(:,end) = solp - relax*bcand(:,end);
%      fcand(end) = .5*socand(:,end)'*Lhs*socand(:,end) + socand(:,end)'*Rhs + ...
%                    mu*(Isum'*abs(Ax*socand(:,end))+Isum'*abs(Ay*socand(:,end)));
%      
%      % Choose the best one
%      [~,icand] = min(fcand);
%      dire = bcand(:,icand);
      
%      dire = Lhs*solp + Rhs + b; % Steepest descent direction
      
%      crit(i) = norm(Bxy*(solp-sol))/norm(Bxy*solp);
%      if i>1 diff(i) = crit(i)-crit(i-1); end
      
%      if i>1 sol = relax*sol + (1-relax)*solp; end % Relaxation step
%      if i>1 sol = solp - relax*dire; end % Steepest descent
      sol = solp - relax*dire;
      
      relaxrec(i)  = relax;
%      quadratic(i) = .5*sol'*Lhs'*Lhs*sol + sol'*Lhs'*Rhs + .5*Rhs'*Rhs;
      quadratic(i) = .5*sol'*Lhs*sol + sol'*Rhs;
      absolute(i)  = mu*(Isum'*abs(Ax*sol)+Isum'*abs(Ay*sol));
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
