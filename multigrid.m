close all;
clear all;

E          = 210000; % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 1;      % N.mm-2 : Loading on the plate
mat        = [0, E, nu];
niter      = 3;
ncase      = 2;
precond    = 2; % 0=no prec; 1=Jacobi; 2=multiscale
init       = 0; % Wether to initialize by a coarse computation

tic

dirichlet = [ 2,1,0 ; 2,2,0 ; 2,3,0 ];
neumann1  = [ 1,3,-fscalar ];
neumann2  = [ 1,1,fscalar ];

[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/multigrid/plate3d_crack0.msh' );
[ nodesc,elementsc,ntoelemc,boundaryc,orderc ] = readmesh3D( 'meshes/multigrid/plate3d_crack1.msh' );

nnodes = size(nodes,1); nnodesc = size(nodesc,1);

Kinter = Krig3D (nodes,elements,mat,order,boundary,[]);
C = Cbound2 ( nodes, dirichlet, boundary ); index = find(diag(C*C'));
indexdex = setdiff(1:3*nnodes,index);
f = [ loading3D(0,nodes,boundary,neumann1), loading3D(0,nodes,boundary,neumann2) ];

Kinterc = Krig3D (nodesc,elementsc,mat,orderc,boundaryc,[]);
Cc = Cbound2 ( nodesc, dirichlet, boundaryc ); indexc = find(diag(Cc*Cc'));
indexdexc = setdiff(1:3*nnodesc,indexc);
fc = [ loading3D(0,nodesc,boundaryc,neumann1), loading3D(0,nodesc,boundaryc,neumann2) ];
uc = zeros(3*nnodesc,ncase);

toc 
tic

%% Direct resolution (for reference)
uref = zeros(3*nnodes,ncase);
uref(indexdex,:) = Kinter(indexdex,indexdex)\f(indexdex,:);

toc
tic

ndofs  = 3*nnodes;
ndofsc = 3*nnodesc;

PJag = 1./diag(Kinter);
PJac = sparse( 1:ndofs, 1:ndofs, PJag, ndofs, ndofs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coarse-preconditionned Conjugate gradient
d     = zeros( ndofs, ncase*(niter+1) );  % ncase directions per step
Ad    = zeros( ndofs, ncase*(niter+1) );
Res   = zeros( ndofs, ncase*(niter+1) );
Zed   = zeros( ndofs, ncase*(niter+1) ); % Yep, we store all that (my name is Leak, Memory Leak)
   
% Equivalent Saad matrices
%unsuralpha   = zeros( ncase*(niter+1) );
%betasuralpha = zeros( ncase*(niter+1) );
%etaeta       = zeros( ncase*(niter+1) );

%fco = passMesh3DF( nodes, elements, nodesc, elementsc, f );
%plotGMSH3D( {fco(:,1),'f2';fc(:,1),'f1'}, elementsc, nodesc, 'output/fco' );
%bug;

if init == 0
   u0 = zeros(3*nnodes,ncase);
else
   u0c(indexdexc,:) = Kinterc(indexdexc,indexdexc)\fc(indexdexc,:);
   u0 = passMesh3Db( nodesc, elementsc, nodes, elements, u0c );
end
u = u0;

Axz = zeros(3*nnodes,ncase);
Axz(indexdex,:) = Kinter(indexdex,indexdex) * u0(indexdex,:);
b = f;

Res(:,[1:ncase]) = b - Axz;

if precond == 0
   Zed(:,[1:ncase]) = Res(:,[1:ncase]);
elseif precond == 1
   Zed(:,[1:ncase]) = PJac*Res(:,[1:ncase]);
else
   fco = passMesh3DF( nodes, elements, nodesc, elementsc, Res(:,[1:ncase]) );
   uco = zeros(ndofsc,ncase);
   uco(indexdexc,:) = Kinterc(indexdexc,indexdexc)\fco(indexdexc,:);
   Zed(:,[1:ncase]) = passMesh3Db( nodesc, elementsc, nodes, elements, uco );
end

d(:,[1:ncase]) = Zed(:,[1:ncase]);
   
residual(1) = norm( Res, 'fro' );
error(1)    = norm(u - uref, 'fro') / norm(uref, 'fro');

%%
V  = zeros(ndofs, ncase*niter);
AV = zeros(ndofs, ncase*niter);
MV = zeros(ndofs, ncase*niter);
H  = zeros(ncase*niter);
   
num = zeros(ncase); % useless, but eta needs initialization #lazy
den = zeros(ncase);

%%
for iter = 1:niter
    multindex   = [ncase*(iter-1)+1:ncase*iter];
    multindexp1 = multindex + ncase;
    multindexm1 = multindex - ncase;

    %% Optimal step
    Ad(indexdex,multindex) = Kinter(indexdex,indexdex) * d(indexdex,multindex);

    denprec = den; numprec = num; % Store those ones
    den = d(:,multindex)'*Ad(:,multindex);

    if rank(den) ~= size(den,1)
       warning('Your test-cases are of lower dimension as their number');
    end

    sqD = den^(1/2);
    d(:,multindex) = d(:,multindex) * inv(sqD);
    Ad(:,multindex) = Ad(:,multindex) * inv(sqD);
    num = Res(:,multindex)'*Zed(:,multindex);
    num = sqD\num; % because of Zed and not d
       
    u = u + d(:,multindex)*num;
    Res(:,multindexp1) = Res(:,multindex) - Ad(:,multindex)*num;
       
    residual(iter+1) = norm( Res(:,multindexp1), 'fro' );
    error(iter+1)    = norm(u-uref,'fro') / norm(uref,'fro');

    if precond == 0
       Zed(:,multindexp1) = Res(:,multindexp1);
    elseif precond == 1
       Zed(:,multindexp1) = PJac*Res(:,multindexp1);
    else
       fco = passMesh3DF( nodes, elements, nodesc, elementsc, Res(:,multindexp1) );
       uco = zeros(ndofsc,ncase);
       uco(indexdexc,:) = Kinterc(indexdexc,indexdexc)\fco(indexdexc,:);
       Zed(:,multindexp1) = passMesh3Db( nodesc, elementsc, nodes, elements, uco );
    end

%    % Ritz variables (Saad++)
%    unsuralpha( multindex, multindex ) = (sqD*num)^(-1/2)*den*(sqD*num)^(-1/2);
%       
%    if iter > 1
%       betasuralpha( multindexm1, multindexm1 ) = ...
%                         (sqD*num)^(-1/2) * betaij'*betaij * (sqD*num)^(-1/2);
%                         % use betaij from the previous iteration
%                            
%       etaeta( multindex, multindex ) = ...
%               (denprec^(1/2)*numprec)^(-1/2) * denprec * ...
%               inv( denprec^(1/2)*numprec ) * (sqD*num)^(1/2);
%    end
       
%    % Reorthogonalize the residual (as we use it next), in sense of M
%    for jter=1:iter-1
%        multjndex = [ncase*(jter-1)+1:ncase*jter];
%        betac = (Res(:,multjndex)'*Zed(:,multjndex)) \...
%                (Res(:,multjndex)'*Zed(:,multindexp1)) ;
%   
%        Zed(:,multindexp1) = Zed(:,multindexp1) - Zed(:,multjndex) * betac;
%        Res(:,multindexp1) = Res(:,multindexp1) - Res(:,multjndex) * betac;
%    end
   
    %% Orthogonalization
    d(:,multindexp1) = Zed(:,multindexp1);
    for jter=iter:iter % Only the last one
        multjndex = [ncase*(jter-1)+1:ncase*jter];
        betaij = ( Ad(:,multjndex)'*d(:,multjndex) ) \ ...
            ( Ad(:,multjndex)'*d(:,multindexp1) );
   
        d(:,multindexp1) = d(:,multindexp1) - d(:,multjndex) * betaij;
    end
   
%    %% The Ritz elements
%    V(:,multindex) = (-1)^(iter-1) * Zed(:,multindex) * ...
%      ((Res(:,multindex)'*Zed(:,multindex))^(-1/2)) ;
%                          
%    delta  = unsuralpha( multindex, multindex ) ;
%    if iter > 1
%       delta = delta + betasuralpha( multindexm1, multindexm1 );
%    end
%   
%    eta = etaeta( multindex, multindex ); % what a stupid variable name
%       
%    if iter > 1
%       H( multindex , [multindexm1,multindex] ) = [eta', delta];
%       H( multindexm1 , multindex ) = eta;
%    else
%       H( multindex , multindex ) = delta;
%    end
end

% Compute eigenelems of the Hessenberg :
%[Q,Theta1] = eig(H);

toc

figure;
hold on;
loglog(residual/residual(1),'Color','blue');
loglog(error,'Color','red');
legend('residual','error')


toplot = { uref(:,1), 'U1_REF' ; uref(:,2), 'U2_REF' ;...
           u(:,1), 'U1' ; u(:,2), 'U2' };
plotGMSH3D( toplot, elements, nodes, 'output/solution' );
