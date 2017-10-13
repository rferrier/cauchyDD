% 10/10/2017
% Algo Steklov-Poincaré primal avec Gradient Conjugué, cas 3D sphérique

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 15;
precond = 1;      % 1 : Use a dual precond, 2 : use H1/2 precond, 3 : use gradient precond
br      = 0.;    % noise
brt     = 0;      % "translation" noise
ntrunc  = 0;     % Nb of Ritz modes (+1)
%inhomog = 0;      % inhomogeneous medium
disc_me = 0;      % Use only a limited number of points as measure
disc_id = 1;      % Identify only a limited number of points
R2      = 15;     % Exterior radius (gemoetrical data)
R1      = 5;      % Interior radius

mespoint = ones(8,8); % Gives the measure points in case disc_me==1
% mespoint : raw = theta, column = phi
%thetames = [1,0,1,0,1,0,1,0]; phimes   = [1,0,1,0,1,0,1,0];
%thetames = [1,1,1,1,1,0,0,0]; phimes   = [1,1,1,1,1,0,0,0];
%mespoint = thetames'*phimes;

idepoint = ones(8,8); % Gives the measure points in case disc_me==1
% mespoint : raw = theta, column = phi
%thetaid = [1,0,1,0,1,0,1,0]; phiid   = [1,0,1,0,1,0,1,0];
thetaid = [1,1,1,1,1,0,0,0]; phiid   = [1,1,1,1,1,0,0,0];
idepoint = thetaid'*phiid;

% methods : 1-> SPP
methods = [1];

mat = [0, E, nu];

dirichlet = [ 0,1,0 ; 0,2,0 ; 0,3,0 ; 0,4,0 ; 0,5,0 ; 0,6,0 ];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/sphere/sphere.msh' );
nnodes = size(nodes,1);

% Then, build the stiffness matrix :
[Kinter,~,~] = Krig3D (nodes,elements,mat,order,boundary,[]);
[ C, node2c, c2node, nbloq ] = Cbound ( nodes, dirichlet, boundary );
K = [Kinter,C;C',zeros(nbloq)];

% The right hand side :
f = pressureLoad3D( nbloq, nodes, boundary, fscalar*[10,0,0,0,1], 2 );

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:3*nnodes,1);
fref = Kinter*uref;
%lagr = uin(3*nnodes+1:end,1);

sigma = stress3D(uref,mat,nodes,elements,order,1,ntoelem);
plotGMSH3D({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'output/reference');

%% Extract the boundary 2
bounindex = find( boundary(:,1)==2 );
patch2    = boundary(bounindex,2:end);
bounindex = find( boundary(:,1)==1 );
patch1    = boundary(bounindex,2:end);

% find the index of the upper nodes
Ckk = Cbound ( nodes, [1,1,0;1,2,0;1,3,0], boundary );
indexC= sum(Ckk'); index1 = find(indexC);
Ckk = Cbound ( nodes, [2,1,0;2,2,0;2,3,0], boundary );
indexC= sum(Ckk'); index2 = find(indexC);

indexa = index2; % Other name

noises = load('./noises/noiseSphere.mat'); % Particular noise vector
noise  = noises.noise;
%noise = randn(size(index1,2),1);
urefb = uref;
urefb(index1) = ( 1 + brt + br*noise ) .* uref(index1);
% H1 norm
Kreg = H1 (nodes, patch2, order); Kreg = Kreg(indexa,indexa);

if disc_me == 1 % Get the discrete measure points
   theta = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4];
   phi = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4];
   X = R2*sin(theta')*cos(phi);
   Y = R2*sin(theta')*sin(phi);
   Z = R2*cos(theta')*ones(1,8);
   nmin = zeros(8,8);
   for i=1:8
      for j=1:8 % For each point, find the closest one in the mesh
         if mespoint(i,j) == 1
            x = X(i,j); y = Y(i,j); z = Z(i,j);
            nodiff = nodes - [x,y,z];
            normdiff = nodiff(:,1).^2+nodiff(:,2).^2+nodiff(:,3).^2;
            [~,nmin(i,j)] = min(normdiff);
         end
      end
   end
   nmin1 = nmin(:); nmin1 = nmin1(find(nmin1)); nmin1 = unique(nmin1);
   indexmeas = [ 3*nmin1(:)-2 ; 3*nmin1(:)-1 ; 3*nmin1(:) ];
   [ Cmeas,~,~, ntotmeas ] = Cnodes ( nodes, indexmeas );
end

if disc_id == 1 % Get the discrete identification points
   theta = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4];
   phi = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4];
   X = R1*sin(theta')*cos(phi);
   Y = R1*sin(theta')*sin(phi);
   Z = R1*cos(theta')*ones(1,8);
   nmin = zeros(8,8);
   for i=1:8
      for j=1:8 % For each point, find the closest one in the mesh
         if idepoint(i,j) == 1
            x = X(i,j); y = Y(i,j); z = Z(i,j);
            nodiff = nodes - [x,y,z];
            normdiff = nodiff(:,1).^2+nodiff(:,2).^2+nodiff(:,3).^2;
            [~,nmin(i,j)] = min(normdiff);
         end
      end
   end
   nmin1 = nmin(:); nmin1 = nmin1(find(nmin1)); nmin1 = unique(nmin1);
   indexid = [ 3*nmin1(:)-2 ; 3*nmin1(:)-1 ; 3*nmin1(:) ];
   [ Cid,~,~, ntotid ] = Cnodes ( nodes, indexid );
end
% Compute norm of displacement
unoref = sqrt(uref(3:3:end).^2+uref(2:3:end-1).^2+uref(1:3:end-2).^2);
urefa = uref(indexa); uref1 = uref(index1);
unorefa = sqrt(urefa(3:3:end).^2+urefa(2:3:end-1).^2+urefa(1:3:end-2).^2);
unoref1 = sqrt(uref1(3:3:end).^2+uref1(2:3:end-1).^2+uref1(1:3:end-2).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

if disc_me == 1
   % First problem
   dirichlet1 = [2,1,0;2,2,0;2,3,0];
   [ C1, ~, ~, nbloq1 ] = Cbound ( nodes, dirichlet1, boundary );
   C1 = [C1,Cmeas]; nbloq1 = nbloq1+ntotmeas;
   K1 = [Kinter,C1;C1',zeros(nbloq1)];
   C12 = Cbound ( nodes, [2,1,0;2,2,0;2,3,0], boundary );
   C11 = Cmeas;
   
   % Second problem
   dirichlet2 = [2,1,0;2,2,0;2,3,0];
   [ C2, ~, ~, nbloq2 ] = Cbound ( nodes, dirichlet2, boundary );
   K2 = [Kinter,C2;C2',zeros(nbloq2)];
   C22 = Cbound ( nodes, [2,1,0;2,2,0;2,3,0], boundary ); % OK, I admit it's the same as C12
   C21 = Cmeas;
   
   %% Dual problems
   % First problem
   dirichlet1d = [];
   C1d = Cmeas; nbloq1d = ntotmeas;
   K1d = [Kinter,C1d;C1d',zeros(nbloq1d)];
   
   % Second problem
   dirichlet2d = [];
   C2d = zeros(3*nnodes,0); K2d = Kinter; nbloq2d = 0;

elseif disc_id == 1
   % First problem
   dirichlet1 = [1,1,0;1,2,0;1,3,0];
   [ C1, ~, ~, nbloq1 ] = Cbound ( nodes, dirichlet1, boundary );
   C1 = [C1,Cid]; nbloq1 = nbloq1 + ntotid;
   K1 = [Kinter,C1;C1',zeros(nbloq1)];
   C12 = Cid;
   C11 = Cbound ( nodes, [1,1,0;1,2,0;1,3,0], boundary );
   
   % Second problem
   dirichlet2 = [];
   C2 = Cid; nbloq2 = ntotid;
   K2 = [Kinter,C2;C2',zeros(nbloq2)];
   C22 = Cid; % OK, I admit it's the same as C12
   C21 = Cbound ( nodes, [1,1,0;1,2,0;1,3,0], boundary );
   
   %% Dual problems
   % First problem
   dirichlet1d = [1,1,0;1,2,0;1,3,0];
   [ C1d, ~, ~, nbloq1d ] = Cbound ( nodes, dirichlet1d, boundary );
   K1d = [Kinter,C1d;C1d',zeros(nbloq1d)];
   
   % Second problem
   dirichlet2d = [];
   C2d = zeros(3*nnodes,0); K2d = Kinter; nbloq2d = 0;
   
else
   % First problem
   dirichlet1 = [2,1,0;2,2,0;2,3,0;
                 1,1,0;1,2,0;1,3,0];
   
   [ C1, ~, ~, nbloq1 ] = Cbound ( nodes, dirichlet1, boundary );
   K1 = [Kinter,C1;C1',zeros(nbloq1)];
   C12 = Cbound ( nodes, [2,1,0;2,2,0;2,3,0], boundary );
   C11 = Cbound ( nodes, [1,1,0;1,2,0;1,3,0], boundary );
   
   % Second problem
   dirichlet2 = [2,1,0;2,2,0;2,3,0];
   [ C2, ~, ~, nbloq2 ] = Cbound ( nodes, dirichlet2, boundary );
   K2 = [Kinter,C2;C2',zeros(nbloq2)];
   C22 = Cbound ( nodes, [2,1,0;2,2,0;2,3,0], boundary ); % OK, I admit it's the same as C12
   C21 = Cbound ( nodes, [1,1,0;1,2,0;1,3,0], boundary );
   
   %% Dual problems
   % First problem
   dirichlet1d = [1,1,0;1,2,0;1,3,0];
   [ C1d, ~, ~, nbloq1d ] = Cbound ( nodes, dirichlet1d, boundary );
   K1d = [Kinter,C1d;C1d',zeros(nbloq1d)];
   
   % Second problem
   dirichlet2d = [];
   C2d = zeros(3*nnodes,0); K2d = Kinter; nbloq2d = 0;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
if find(methods==1)

   % First, plot the reference
   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',unoref,'FaceColor','interp');
   colorbar;
   
   figure;
   ret = patch('Faces',patch1,'Vertices',nodes, ... 
               'FaceVertexCData',unoref,'FaceColor','interp');
   colorbar;
   
   error    = zeros(niter+1,1);
   residual = zeros(niter+1,1);
   regulari = zeros(niter+1,1);
   
   Itere  = zeros( 3*nnodes, 1 );
   d      = zeros( 3*nnodes, niter+1 );
   Ad     = zeros( 3*nnodes, niter+1 );
   Res    = zeros( 3*nnodes, niter+1 );
   Zed    = zeros( 3*nnodes, niter+1 );
   alpha  = zeros( niter+1, 1 );
   beta   = zeros( niter+1, 1 );
   alpha2 = zeros( niter+1, 1 );

   % Pre-computation
   C1tC12C12t = C1'*C12*C12';
   C12C12tC1  = C12*C12'*C1;
   C2tC22C22t = C2'*C22*C22';
   C22C22tC2  = C22*C22'*C2;
   C12C12t    = C12*C12';
   
   %% Perform A x0 :
   % Solve 1
   f1 = [ zeros(3*nnodes,1) ; C1tC12C12t*Itere ]; % u=0 outside of 2
   uin1 = K1\f1;
   lagr1 = uin1(3*nnodes+1:end,1);
   lamb1 = -C12C12tC1*lagr1;    % restriction of lagr to 2
   % Solve 2
   f2 = [ zeros(3*nnodes,1) ; C2tC22C22t*Itere ];
   uin2 = K2\f2;
   lagr2 = uin2(3*nnodes+1:end,1);
   lamb2 = -C22C22tC2*lagr2;
   %
   Axz = lamb1-lamb2;
   %%%%
   %% Compute Rhs :
   % Solve 1
   f1 = [ zeros(3*nnodes,1) ; C1'*C11*C11'*urefb ];
   uin1 = K1\f1;
   lagr1 = uin1(3*nnodes+1:end,1);
   lamb1 = -C12C12tC1*lagr1;
   % Solve 2 (zero 'cause no loading)
   f2 = loading3D(nbloq2,nodes,boundary,[]);
   uin2 = K2\f2;
   lagr2 = uin2(3*nnodes+1:end,1);
   lamb2 = -C22C22tC2*lagr2;
   %
   b = lamb2-lamb1;
   %%
   Res(:,1) = b - Axz;
   
   if precond == 1
       % Solve 1
       f1 = [Res(:,1); zeros(nbloq1d,1)];
       uin1 = K1d\f1;
       u1i = uin1(1:3*nnodes,1);
       u1 = C12C12t*u1i;
       Zed(:,1) = u1;
   else
       Zed(:,1) = Res(:,1);
   end
   
   d(:,1) = Zed(:,1);
   
   residual(1) = norm(Res( :,1));
   error(1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
   regulari(1) = sqrt(Itere(indexa)'*Kreg*Itere(indexa));%norm(Itere);
   
   eta      = 0;
   %%
   for iter = 1:niter
       % Solve 1
       f1 = [ zeros(3*nnodes,1) ; C1tC12C12t*d(:,iter) ]; % u=0 outside of 2
       uin1 = K1\f1;
       lagr1 = uin1(3*nnodes+1:end,1);
       lamb1 = -C12C12tC1*lagr1;    % restriction of lagr1 to 2
       % Solve 2
       f2 = [ zeros(3*nnodes,1) ; C2tC22C22t*d(:,iter) ];
       uin2 = K2\f2;
       lagr2 = uin2(3*nnodes+1:end,1);
       lamb2 = -C22C22tC2*lagr2;
       %
        Ad(:,iter) = lamb1-lamb2;
       
       den = (d(:,iter)'*Ad(:,iter));
       d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
       num = Res(:,iter)'*d(:,iter);
       
       Itere         = Itere + d(:,iter)*num;%/den;
       Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
       
       residual(iter+1) = norm(Res(:,iter+1));
       error(iter+1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
       regulari(iter+1) = sqrt(Itere(indexa)'*Kreg*Itere(indexa));%norm(Itere);
       
       if precond == 1
           % Solve 1
           f1 = [Res(:,iter+1); zeros(nbloq1d,1)];
           uin1 = K1d\f1;
           u1i = uin1(1:3*nnodes,1);
           u1 = C12C12t*u1i;
           Zed(:,iter+1) = u1;
       else
           Zed(:,iter+1) = Res(:,iter+1);
       end
       
       % Needed values for the Ritz stuff
       alpha(iter) = num/sqrt(den);
       alpha2(iter) = num;
       beta(iter)  = - Zed(:,iter+1)'*Ad(:,iter)/sqrt(den);
       
       % First Reorthogonalize the residual (as we use it next), in sense of M
       for jter=1:iter-1
           betac = Zed(:,iter+1)'*Res(:,jter) / (Zed(:,jter)'*Res(:,jter));
           Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
       end
       
       %% Orthogonalization
       d(:,iter+1) = Zed(:,iter+1);
       
       for jter=iter:iter % No need to reorthogonalize (see above)
           betaij = ( Zed(:,iter+1)'*Ad(:,jter) );%/...
               %( d(:,jter)'*Ad(:,jter) );
           d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
   
       end
       
       %% Ritz algo : find the Ritz elements
       % Build the matrices
       V(:,iter) = zeros(3*nnodes,1);
       V(:,iter) = (-1)^(iter-1)*Zed(:,iter)/(sqrt(Res(:,iter)'*Zed(:,iter)));
       %norm(Zed(:,iter))
       etap   = eta;
       delta  = 1/alpha(iter);
       if iter > 1
          delta  = delta + beta(iter-1)/alpha(iter-1);
       end
       eta    = sqrt(beta(iter))/alpha(iter);
       
       if iter > 1
          H(iter,[iter-1,iter]) = [etap, delta];
          H(iter-1,iter)        = etap;
       else
          H(iter,iter) = delta;
       end
       
       % Compute eigenelems of the Hessenberg :
       [Q,Theta1] = eig(H);
       theta = diag(Theta1);
       % Sort it
       [theta,Ind] = sort(theta,'descend');
       Q = Q(:,Ind);
       Theta1 = Theta1(Ind,Ind);
       Y = V*Q;
   end
   
   % hack to avoid the warnings (theta1 is diagonal
   Itheta1 = diag(1./diag(Theta1));
   
   % Compute the solution
   chi = Itheta1*Y'*b;
   if ntrunc > 0
      chi(ntrunc:end) = 0;
   end
   ItereR = Y*chi;
   
   %%for i=1:niter
   %   plot(V(2*b2node2-1,20:25))
   %   figure;
   %%end
   
   regS = zeros(niter,1);
   resS = zeros(niter,1);
   %% Build the L-curve regul, ntrunc
   for i = 1:iter+1
      chiS   = Itheta1*Y'*b; chiS(i:end) = 0;
      ItereS = Y*chiS;
      
%      % Solve 1
%      f1 = [ zeros(3*nnodes,1) ; C1tC12C12t*ItereS ]; % u=0 outside of 2
%      uin1 = K1\f1;
%      lagr1 = uin1(3*nnodes+1:end,1);
%      lamb1 = -C12C12tC1*lagr1;    % restriction of lagr to 2
%      % Solve 2
%      f2 = [ zeros(3*nnodes,1) ; C2tC22C22t*ItereS ];
%      uin2 = K2\f2;
%      lagr2 = uin2(3*nnodes+1:end,1);
%      lamb2 = -C22C22tC2*lagr2;
%      %
%      AI = lamb1-lamb2;
%      
%      ResS = AI-b;
%      resS(i) = norm(ResS);   
%      regS(i) = sqrt(ItereS(indexa)'*Kreg*ItereS(indexa));%norm(ItereS);
      errS(i) = norm(ItereS(indexa)-uref(indexa))/norm(uref(indexa));
   end
   
   % Residual in the diagonal basis :
   %resD = zeros(iter,1);  regD = zeros(iter,1);  bt = Y'*b;
   %for i=1:iter
   %   resD(i) = sqrt( sum( bt(i:end).^2) );
   %   regD(i) = sqrt( sum((theta(1:i) .* chit(1:i)).^2) );
   %end
   traD = 1 - sum(theta(1:ntrunc-1))/sum(theta);
   regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
   for i = 1:iter+1
      chiD   = Itheta1*Y'*b; chiD(i:end) = 0;
      ItereD = Y*chiD;
      %AI = Y*Theta1*Y'*ItereD;
      %ResD = AI-b;
      %resD(i) = sqrt(norm(ResD));  %Problem with the preconditionner
      resD(i) = sqrt( sum( bt(i:end).^2) );  
      regD(i) = sqrt(ItereD(indexa)'*Kreg*ItereD(indexa));%norm(ItereD);
   end
   
   [indm2,pol] = findPicard2(log10(abs(chiD)), 5, 1);
   n = size(pol,1);
   t = 1:.05:niter; tt = zeros(n,20*(niter-1)+1);
   for j=1:n
      tt(j,:) = t.^(n-j);
   end
   px = pol'*tt;

   figure;
   hold on;
   plot(log10(theta),'Color','blue')
   plot(log10(abs(Y'*b)),'Color','red')
   plot(log10(abs(chiD)),'Color','black')
   plot(t,px,'Color','cyan')
   legend('Ritz Values','RHS values','solution coefficients', ...
          'polynomial interpolation')

   %
%   errorI = (Itere(3:3:end)-uref(3:3:end))/norm(uref(3:3:end));
%   nerrorI = norm(Itere(indexa)-uref(indexa))/norm(uref(indexa));
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',errorI,'FaceColor','interp');
%   colorbar;
%   %
%   errorR = (ItereR(3:3:end)-uref(3:3:end))/norm(uref(3:3:end));
%   nerrorR = norm(ItereR(indexa)-uref(indexa))/norm(uref(indexa));
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',errorR,'FaceColor','interp');
%   colorbar;
   %
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',Itere(3:3:end),'FaceColor','interp');
%   colorbar;
%   %

   % Compute norm of solution displacement
   Itereno = sqrt(ItereR(3:3:end).^2+ItereR(2:3:end-1).^2+ItereR(1:3:end-2).^2);

   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',Itereno,'FaceColor','interp');
   colorbar;
   
   relerror = abs(Itereno-unoref)/max(abs(unorefa));
   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',relerror,'FaceColor','interp');
   colorbar;
   
   figure;
   hold on
   plot(error(2:end),'Color','blue')
%   plot(residual(2:end),'Color','red')
   plot(errS(2:end),'Color','green')
   legend('error','Ritz error')
%   legend('error')
   
%   %L-curve :
%   figure;
%   hold on;
%   set(gca, 'fontsize', 20);
%   loglog(residual(2:iter+1),regulari(2:iter+1),'Color','red','-*','linewidth',3);
%   loglog(resS(2:iter+1),regS(2:iter+1),'-+','linewidth',3);
%   legend('L-curve','RL-curve')
%   xlabel('Residual')
%   ylabel('H1 norm')
%   %legend('RL-curve')
%   %findCorner (residual(2:iter+1), regulari(2:iter+1), 3)
%   %findCorner (resS(2:iter+1), regS(2:iter+1), 3)
   
   %%hold on;
   %%loglog(resS(2:iter+1),regS(2:iter+1),'Color','red','-*');
   %loglog(resD(2:iter),regD(2:iter),'-+');
   %%legend('RL-curve (natural basis)','RL-curve (diagonal basis)')
   %figure
   %findCorner (resD(2:iter), regD(2:iter), 3)
   %
   %
   %hold on;
   %plot(fsol(2*b2node2-1),'Color','red')
   %plot(fsolR(2*b2node2-1),'Color','blue')
   %plot(f(2*b2node2-1),'Color','green')
   %legend('brutal solution','filtred solution','reference')
   %
   %total_error = norm(usol-uref)/norm(uref);
   %% Compute stress :
   %sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
   %plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
   
   %% Compute the final solution
   ffinal = [ zeros(3*nnodes,1) ; C2tC22C22t*ItereR ];
   uinfinal = K2\ffinal; ufinal = uinfinal(1:3*nnodes);
   
   ufinalno = sqrt(ufinal(3:3:end).^2+ufinal(2:3:end-1).^2+ufinal(1:3:end-2).^2);
   
   figure;
   ret = patch('Faces',patch1,'Vertices',nodes, ... 
               'FaceVertexCData',ufinalno,'FaceColor','interp');
   colorbar;
   
   relerrore = abs(ufinalno-unoref)/max(abs(unoref1));
   figure;
   ret = patch('Faces',patch1,'Vertices',nodes, ... 
               'FaceVertexCData',relerrore,'FaceColor','interp');
   colorbar;
   errorext = norm(ufinal(index1)-uref(index1))/norm(uref(index1));
end