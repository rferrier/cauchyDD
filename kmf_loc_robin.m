%13/09/2016
%Algo KMF-Robin avec localisation

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E        = 70000;  % MPa : Young modulus
nu       = 0.3;    % Poisson ratio
fscalar  = 1;      % N.mm-1 : Loading on the plate
niter    = 5;
br       = 0.5;      % noise
relax    = 0;      % Use a relaxation paramter
nlociter = 1;      % Nb of localization iterations
max_erc  = 2;      % Erc criterion for localization
k0       = E;      % Basic Robin multiplicator
Lc       = 0;      % Correlation length for the white noise

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
[ node2b2, b2node2 ] = mapBound( 2, boundary, nnodes );
index    = 2*b2node3-1;
index    = index(size(index):-1:1);
index2   = 2*b2node2-1;

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

%brbru = randn(2*nnodes,1);
%brtro = correfilter (nodes, Lc, 'rectangle', brbru);
%hold on;
%plot( brbru(index2), 'Color', 'blue' );
%plot( brtro(index2), 'Color', 'red' );

urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
%urefb = correfilter (nodes, Lc, 'rectangle', urefb);
%urefb = ( 1 + br*brtro ) .* uref;
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

error1   = zeros(niter*nlociter,1);
error2   = zeros(niter*nlociter,1);
residual = zeros(niter*nlociter,1);
regulari = zeros(niter*nlociter,1);
serror1   = zeros(niter,1); % Error for sigma
serror2   = zeros(niter,1);
sresidual = zeros(niter,1);
ferror1   = zeros(niter,1); % Error for reaction force
ferror2   = zeros(niter,1);
fresidual = zeros(niter,1);
Erdc      = ones(niter*nlociter,1);
Erdc1     = ones(niter*nlociter,1);
Erdc2     = ones(niter*nlociter,1);

% init :
u1 = uref-uref;
u2 = u1; fri= u1; fr2= u1; v  = u1;
theta = ones(niter+1,1); % First relaxation parameter
Ko   = k0*ones(size(boundary,1),1);% First Robin parameter

% DN problem
dirichlet1 = [4,1,0;4,2,0;
             3,1,0;3,2,0];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);

% ND problem
dirichlet2 = [4,1,0;4,2,0];
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);

for Bigi = 1:nlociter
  
  for iter = 1:niter
      % Solve DN
      u1o = u1;
      f1in = [f1ini ; zeros(nbloq1,1)];
      fdir = dirichletRhs2(v, 3, c2node1, boundary, nnodes );
      f1 = f1in + fdir;
  
      uin1 = K1\f1;
      u1 = uin1(1:2*nnodes,1);
      lagr1 = uin1(2*nnodes+1:end,1);
      frio = fri; % Store fri for the residual computation
      fri = Kinter*u1;
      
  %     sigma1 = stress(u1,E,nu,nodes,elements,order,1,ntoelem);
      error1(niter*(Bigi-1)+iter) = sqrt((u1-uref)'*(u1-uref)/(uref'*uref));
      ferror1(niter*(Bigi-1)+iter) = sqrt((fri-fref)'*(fri-fref)/(fref'*fref));
      
      % Solve ND
      u2o = u2;
      fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size
      [Krob, fdir1] = robinRHS( nbloq2, nodes, boundary, urefb, Ko, 1 );
      %fdir1 = dirichletRhs2(urefb, 1, c2node2, boundary, nnodes);

      % Assembly the dirichlet RHS (ensure a redundant CL isn't imposed 2
      % times
      f2 = fr + fdir1; %assembleDirichlet( [fdir1,fdir2] );
  
      uin2 = (K2+Krob)\f2;
      u2 = uin2(1:2*nnodes,1);
      fr2o = fr2;
      fr2 = Kinter*u2;
      
      vo = v;
      v = theta(iter)*u2 + (1-theta(iter))*vo;
      
      if relax == 1 && iter > 1
          e1 = u1-u1o;
          e2 = u2-u2o;
          theta(iter+1) = myps(e1,e1-e2,Kinter,boundary,M,nodes) /...
              myps(e1-e2,e1-e2,Kinter,boundary,M,nodes);
      end
      
      error2(niter*(Bigi-1)+iter) = sqrt((u2-uref)'*(u2-uref)/(uref'*uref));
      ferror2(niter*(Bigi-1)+iter) = sqrt((fr2-fref)'*(fr2-fref)/(fref'*fref));
      
      residual(niter*(Bigi-1)+iter) = sqrt((u1-u2)'*(u1-u2)/sqrt((u1'*u1)*(u2'*u2)));
      fresidual(niter*(Bigi-1)+iter) = sqrt((fri-fr2)'*(fri-fr2)/sqrt((fri'*fri)*(fr2'*fr2)));
                   
      regulari(niter*(Bigi-1)+iter) = sqrt(u2'*regul(u2, nodes, boundary, 3));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Computation of the constitutive law error (at each boundary element)
      Erc    = zeros(nbound,1);  % This one stores the error
      Erc1   = zeros(nbound,1);  % real Erc
      Erc2   = zeros(nbound,1);  % real Erc
      nberc  = 0; %Nb of elements that contribute to Erc
      for i=1:nbound
         if boundary(i,1) == 3 % Missing boundary
            no1 = boundary(i,2);
            no2 = boundary(i,3);
            x1 = nodes(no1,1); y1 = nodes(no1,2);
            x2 = nodes(no2,1); y2 = nodes(no2,2);
            %leng = sqrt( (x2-x1)^2 + (y2-y1)^2 );
            map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
            
            Me = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
            deltaf = frio(map) - fri(map);
            deltau = u2(map) - u2o(map);
            %Erc(i) = deltaf'*Me*deltau;
            %nberc = nberc+1;
         end
         if boundary(i,1) == 1 % Redondant data boundary
            no1 = boundary(i,2);
            no2 = boundary(i,3);
            x1 = nodes(no1,1); y1 = nodes(no1,2);
            x2 = nodes(no2,1); y2 = nodes(no2,2);
            %leng = sqrt( (x2-x1)^2 + (y2-y1)^2 );
            map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
            
            Me = [1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1]; % Well, it's regular scalar product...
            deltaf = fr2(map) - fref(map);   %/!\ Noise on f isn't taken into account
            deltau = urefb(map) - u1(map);
            Erc(i) = deltaf'*Me*deltau;
            nberc = nberc+1;

       end
    end
    Erc = abs(Erc); % Why are there < 0 therms ?
    Erdc(niter*(Bigi-1)+iter) = sum(Erc);
  end
  
  % Tune Ko in function of Erc
  Ko = k0*Erdc(niter*(Bigi-1)+iter)./(Erc*nberc);
  
%  figure;
%  plot(Erc);
%  figure;
%  plot(uref-urefb);
%  figure
%  plot(Ko);

end

Erdc = Erdc/Erdc(1);  % In order to normalize the stuff
residual(1) = 1; % tiny hack
figure;
hold on;
set(gca, 'fontsize', 15);
set(gca,'ylim',[-3 0])
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
%plot(log10(residual),'Color','red')
plot(log10(ferror1),'Color','yellow')
plot(log10(ferror2),'Color','cyan')
%plot(log10(fresidual),'Color','magenta')
plot(log10(Erdc),'Color','green')
%legend('error1 (log)','error2 (log)','residual (log)','ferror1 (log)',...
       %'ferror2 (log)','fresidual (log)','Erc(log)')
legend('error1 (log)','error2 (log)','ferror1 (log)',...
       'ferror2 (log)','Erc(log)')

% L-curve
%figure;
%loglog(Erdc,regulari);

% Compute stress :
sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({u2,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'field2');