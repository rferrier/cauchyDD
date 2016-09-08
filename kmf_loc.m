%21/03/2016
%Algo KMF

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E        = 70000;  % MPa : Young modulus
nu       = 0.3;    % Poisson ratio
fscalar  = 1;      % N.mm-1 : Loading on the plate
niter    = 20;
br       = 0.01;      % noise
relax    = 0;      % Use a relaxation paramter
nlociter = 2;      % Nb of localization iterations
max_erc  = 2;     % Erc criterion for localization

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [4,1,0;
             4,2,0];
neumann   = [1,2,fscalar;
             2,1,fscalar;
             3,2,fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);
nbound = size(boundary,1);

% Extract the index of the boundary
[ node2b3, b2node3 ] = mapBound( 3, boundary, nnodes );
index    = 2*b2node3-1;
index    = index(size(index):-1:1);

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);

% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);
urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,1) ) .* lagr;
fref = f( 1:2*nnodes,1 ); % Reaction forces

% Compute stress :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({uref,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'reference');
% patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceAlpha',0);
% axis equal
% figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize redondant boundary
neumann1   = [1,2,fscalar;
              2,1,fscalar];
f1ini = loading(0,nodes,boundary,neumann1);
for i=1:nbound
   if boundary(i) == 2
      boundary(i) = 1;
   end
end

for Bigi = 1:nlociter
  % init :
  u1    = uref-uref;
  u2    = u1;
  fri   = u1;
  v     = u1;
  theta = ones(niter+1,1); % First relaxation parameter
  
  % DN problem
  dirichlet1 = [4,1,0;4,2,0;
                3,1,0;3,2,0];
  [K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
  %[L1,U1] = lu(K1);
  % ND problem
  dirichlet2 = [4,1,0;4,2,0;
                1,1,0;1,2,0;
                2,1,0;2,2,0];
  [K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
  %[L2,U2] = lu(K2);
  
  error1   = zeros(niter,1);
  error2   = zeros(niter,1);
  residual = zeros(niter,1);
  regulari = zeros(niter,1);
  serror1   = zeros(niter,1); % Error for sigma
  serror2   = zeros(niter,1);
  sresidual = zeros(niter,1);
  ferror1   = zeros(niter,1); % Error for reaction force
  ferror2   = zeros(niter,1);
  fresidual = zeros(niter,1);
  
  for iter = 1:niter
      % Solve DN
      u1o = u1;
      f1in = [f1ini ; zeros(nbloq1,1)];
      %fdir = dirichletRhs(u2, 3, C1, boundary);
      fdir = dirichletRhs2(v, 3, c2node1, boundary, nnodes );
      f1 = f1in + fdir;
  
      uin1 = K1\f1;
      u1 = uin1(1:2*nnodes,1);
      lagr1 = uin1(2*nnodes+1:end,1);
      frio = fri; % Store fri for the residual computation
      fri = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );
      
  %     sigma1 = stress(u1,E,nu,nodes,elements,order,1,ntoelem);
      
      error1(iter) = sqrt( myps(u1-uref,u1-uref,Kinter,boundary,M,nodes)/...
          myps(uref,uref,Kinter,boundary,M,nodes) );
      
      % Solve ND
      u2o = u2;
      fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size
      fdir1 = dirichletRhs2(urefb, 1, c2node2, boundary, nnodes);
      %fdir2 = dirichletRhs2(urefb, 2, c2node2, boundary, nnodes);
      % Assembly the dirichlet RHS (ensure a redundant CL isn't imposed 2
      % times
      f2 = fr + fdir1; %assembleDirichlet( [fdir1,fdir2] );
  
      uin2 = K2\f2;
      u2 = uin2(1:2*nnodes,1);
      fr2 = Kinter*u2;
      
      vo = v;
      v = theta(iter)*u2 + (1-theta(iter))*vo;
      
      if relax == 1 && iter > 1
          e1 = u1-u1o;
          e2 = u2-u2o;
          theta(iter+1) = myps(e1,e1-e2,Kinter,boundary,M,nodes) /...
              myps(e1-e2,e1-e2,Kinter,boundary,M,nodes);
      end
  
  %     sigma2 = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
      
      error2(iter) = sqrt( myps(u2-uref,u2-uref,Kinter,boundary,M,nodes)/...
          myps(uref,uref,Kinter,boundary,M,nodes) );
      
      residual(iter) = sqrt( myps(u1-u2,u1-u2,Kinter,boundary,M,nodes)/...
                       sqrt(myps(u1,u1,Kinter,boundary,M,nodes)*...
                       myps(u2,u2,Kinter,boundary,M,nodes)) );
                   
      regulari(iter) = sqrt(u2'*regul(u2, nodes, boundary, 3));
  end
  
  residual(1) = 1; % tiny hack
  
  figure;
  hold on;
  set(gca, 'fontsize', 15);
  set(gca,'ylim',[-3 0])
  plot(log10(error1),'Color','black')
  plot(log10(error2),'Color','blue')
  plot(log10(residual),'Color','red')
  legend('error1 (log)','error2 (log)','residual (log)')
  
  %% Computation of the constitutive law error (at each boundary element)
  Erc    = zeros(nbound,1);  % This one stores the error
  for i=1:nbound
     if boundary(i,1) == 3 % Missing boundary
        no1 = boundary(i,2);
        no2 = boundary(i,3);
        x1 = nodes(no1,1); y1 = nodes(no1,2);
        x2 = nodes(no2,1); y2 = nodes(no2,2);
        leng = sqrt( (x2-x1)^2 + (y2-y1)^2 );
        
        Me = leng/3 * [1,0,.5,0;0,1,0,.5;.5,0,1,0;0,.5,0,1];
        deltaf = fri([2*no1-1, 2*no1, 2*no2-1, 2*no2])...
                 - fr2([2*no1-1, 2*no1, 2*no2-1, 2*no2]);
        deltau = u1([2*no1-1, 2*no1, 2*no2-1, 2*no2])...
                 - u2([2*no1-1, 2*no1, 2*no2-1, 2*no2]);
        Erc(i) = deltaf'*Me*deltau;
     end
     if boundary(i,1) == 1 % Redondant data boundary
        no1 = boundary(i,2);
        no2 = boundary(i,3);
        x1 = nodes(no1,1); y1 = nodes(no1,2);
        x2 = nodes(no2,1); y2 = nodes(no2,2);
        leng = sqrt( (x2-x1)^2 + (y2-y1)^2 );
        
        Me = leng/3 * [1,0,.5,0;0,1,0,.5;.5,0,1,0;0,.5,0,1];
        deltaf = fr2([2*no1-1, 2*no1, 2*no2-1, 2*no2])...
                 - fref([2*no1-1, 2*no1, 2*no2-1, 2*no2]);   %/!\ Noise on f isn't taken into account
        deltau = urefb([2*no1-1, 2*no1, 2*no2-1, 2*no2])...
                 - u1([2*no1-1, 2*no1, 2*no2-1, 2*no2]);
        Erc(i) = deltaf'*Me*deltau;
     end
  end
  Erc = abs(Erc); % Why are there < 0 therms ?
  
  figure
  plot(Erc)
  
  % Who will be localized ?
  Erlim = mean(Erc) * max_erc;
  for i=1:nbound
     if (boundary(i,1) == 3 || boundary(i,1) == 1)  % Other ones are safe (for now)
        if Erc(i) > Erlim;
           boundary(i,1) = 3;  % this one is suspect
        else
           boundary(i,1) = 1; % this one is clean
        end
     end
  end
  
  % Replace urefb and f1ini TODO : be sure there are the good ones
  urefb = u2;
  f1ini = fr2;
  
  % L-curve
  %loglog(residual,regulari);
  %figure
end
% Compute stress :
sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({u2,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'field2');