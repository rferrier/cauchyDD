function [u,sig,p,lambda] = Elastoplast( mat,K,f,T,nodes,elem,order,eps,nmax,varargin )
 % This function computes the solution of an elastoplastic problem
 % this was shamelessly copied from "Analyse des solides déformables par la 
 % méthode des éléments finis" by Marc Bonnet and Attilio Frangi.
 %
 % input : mat      : Material properties
 %         K        : Tangent elastic matrix
 %         f        : external forces
 %         T        : time steps
 %         nodes    : List of nodes
 %         elem     : List of elements of the mesh
 %         order    : Order of the elements
 %         eps      : convergence criterion
 %         niter    : maximal number of iterations
 %         order    : order of the elements
 %         varargin : if 1 : plane deformations else plane constrains
 %         varargin : if 1 : display convergence info
 
 % output : u       : solution for every t

 mo = 0;
 if numel(varargin)>0
     mo = cell2mat(varargin(1));
 end
 
 if mo == 0
    warning('Plasticity in plane constraints not implemented yet')
 end
 
 talk = 1;
 if numel(varargin)>1
     talk = cell2mat(varargin(2));
 end

 % Useful stuff
 nnodes = size(nodes,1);
 nelem  = size(elem,1);
 nt     = max( size(T) );
 sizf   = size(f,1);
 nbloq  = sizf - 2*nnodes;
 E = mat(2); nu = mat(3); Slim = mat(4); alpha = mat(5); H = mat(6);
 mu = E/(2*(1+nu)); kappa = E/(3*(1-2*nu));
 sqrt32 = sqrt(3/2); % 'cause time is money
 C = K( 2*nnodes+1:end, 1:2*nnodes); % Dirichlet BCs
 Pro = C'*inv(C*C')*C; % Projector on the constrianed dofs
 Deviator = [ 2/3 -1/3 0 ; -1/3 2/3 0 ; 0 0 1/2 ];% or 1 ?
 
 % Build the elastic behaviour matrix
 if mo == 0
    S = 1/E*[1,-nu,0 ; -nu,1,0 ; 0,0,2*(1+nu)];
    Sm1 = inv(S);
 else
    Sm1 = [4*mu/3+kappa, kappa-2*mu/3, 0;...
           kappa-2*mu/3, 4*mu/3+kappa, 0;...
           0, 0, mu];
 end
 
 % Sanity checks
 if min( size(T,1) ) > 1
    warning('Time step vector is a matrix')
 end
 if size(f,2) ~= nt
    warning('Loading size does not match the number of time steps')
 end
 % Notation : we add f(0) that is 0
 % and      : u(1) corresponds to u0, u(nt+1)  corresponds to the last increment
 f = [zeros(sizf,1) , f];
 % Project f to keep only loads that will actually work with displacement
 f(1:2*nnodes,:) = f(1:2*nnodes,:) - Pro*f(1:2*nnodes,:);
 
 % Internal forces don't act on Dirichlet's bcs : we need incremental Dirichlet conditions
 if nbloq > 0 % How is the contrary possible ?
    f( 2*nnodes+1:end, 2:end ) = f( 2*nnodes+1:end, 2:end ) - f( 2*nnodes+1:end, 1:end-1 );
 end
 
 %% initialization
 u     = zeros( sizf, nt+1 );
 if order == 1
    sigeq = zeros( nelem, nt+1 );
    p     = zeros( nelem, nt+1 );
    X     = zeros( 4*nelem, nt+1 );
    epsp  = zeros( 3*nelem, nt+1 );
    sig   = zeros( 3*nelem, nt+1 );
    sigd  = zeros( 4*nelem, nt+1 ); % Deviator of the previour one (not obvious if plane deformations)
 end
 
 fint   = zeros( sizf, 1 );
 lambda = zeros( 2*nnodes, nt+1 );
 Ktan   = K;
 
 %% Loop over time steps
 for i = 1:nt
    Deltau = zeros( sizf, 1 );
    lambda(:,i+1) = lambda(:,i); % Initialize the reaction forces
    fi     = f(:,i+1) - fint - [ lambda(:,i+1) ; zeros(nbloq,1) ]; % Rhs for this time step
    %fi(1:2*nnodes) = fi(1:2*nnodes) - Pro*fi(1:2*nnodes); % Suppress on bcs
%    resref = norm( fi(1:2*nnodes) )
    stay = 1;
    niter = 1;
    
    sig( : , i+1 )  = sig( : , i ); % Initialize constraint
    sigd( : , i+1 ) = sigd( : , i );
    
    % Loop over Newton steps
    while stay
       % Linear system and Reaction forces
       deltau = Ktan\fi; %norm(deltau)
%       if niter == 1 % The first iteration carries the Dirichlet BC
          lambda(:,i+1) = lambda(:,i+1) - Pro*Ktan(1:2*nnodes,1:2*nnodes)*deltau(1:2*nnodes);
%       end

        if niter == 1 % Recompute resref
           resref = norm( f(:,i+1) - fint - [ lambda(:,i+1) ; zeros(nbloq,1) ] );
        end

%       fmon = Ktan(1:2*nnodes,1:2*nnodes)*deltau(1:2*nnodes);
%       lam2 = lambda(2)
%       lambda(8)
%       
       Deltau = Deltau + deltau;
       Deltap = zeros(size(p,1),1); % Deltap is recomputed at each Newton step
       DeltaX = zeros(size(X,1),1);
       
       % Loop over the elements to build the tangent matrix and the internal forces
       Ktan(1:2*nnodes,1:2*nnodes) = sparse( 2*nnodes, 2*nnodes );
       fint                        = zeros( sizf, 1 );
       for j = 1:nelem
          if order == 1
             n1 = elem(j,1); n2 = elem(j,2); n3 = elem(j,3);
             x1 = nodes(n1,1); y1 = nodes(n1,2);
             x2 = nodes(n2,1); y2 = nodes(n2,2);
             x3 = nodes(n3,1); y3 = nodes(n3,2);
             map = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n3-1, 2*n3];
             
             Deltaue = Deltau( map ); uei = u(map,i);
             
             % There is only 1 Gauss point per element
             S = abs( .5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)) );
             Be=[y2-y3,0,y3-y1,0,y1-y2,0;
                 0,x3-x2,0,x1-x3,0,x2-x1;
                 x3-x2,y2-y3,x1-x3,y3-y1,x2-x1,y1-y2]/(2*S);
             Deltaeps = Be * Deltaue;
             
             %% Radial algo
             Deltaepsd = [ 2/3*Deltaeps(1) - 1/3*Deltaeps(2)
                           2/3*Deltaeps(2) - 1/3*Deltaeps(1)
                           -1/3*Deltaeps(1) - 1/3*Deltaeps(2) % deviator(e33)is not 0
                           Deltaeps(3)/2 ]; % Deviator of Deltaeps (plane defo, eps33=0)
             sigde = sigd( [4*j-3,4*j-2,4*j-1,4*j], i );
             sigdelas = sigde + 2*mu*Deltaepsd; 
             sigeqelas = sqrt32 * sqrt(sigdelas(1)^2 + sigdelas(2)^2 + sigdelas(3)^2 + 2*sigdelas(4)^2);
%             if j==1 sigeqelas
%             Slim + alpha*p(j,i) 
%             Slim + alpha*p(j,i) + alpha*Deltap(j) end
             crit = sigeqelas - (Slim + alpha*p(j,i));
             epsi = Be*uei; epsipp = epsi + Deltaeps;
             if crit <= 0 % No plasticity at this Gauss point
                sig( [3*j-2,3*j-1,3*j], i+1 ) = ...
                      kappa*(epsipp(1)+epsipp(2)) * [1;1;0]... % Hydrostatic part
                      + sigdelas([1,2,4]);
                sigd( [4*j-3,4*j-2,4*j-1,4*j], i+1 ) = sigdelas;
                epsp( [3*j-2,3*j-1,3*j], i+1 ) = epsp( [3*j-2,3*j-1,3*j], i );
                
                Atan = Sm1; % We are in elasticity
             else % Plasticity
                Deltap(j) = ( sigeqelas - Slim - alpha*p(j,i) ) / (alpha+3*mu);
                beta  = 3*mu*Deltap(j)/sigeqelas;
                gamma = 3*mu/(3*mu+alpha);
                sig( [3*j-2,3*j-1,3*j], i+1 ) = ...
                      kappa*(epsipp(1)+epsipp(2)) * [1;1;0]...
                      + (1-beta)*sigdelas([1,2,4]) ;% (1-beta)/2*sigdelas(4) ] ;
                sigd( [4*j-3,4*j-2,4*j-1,4*j], i+1 ) = (1-beta)*sigdelas;
                epsp( [3*j-2,3*j-1,3*j], i+1 ) = epsp( [3*j-2,3*j-1,3*j], i ) + ...
                      beta/(2*mu)*sigdelas([1,2,4]);
                      
%                Atan = Sm1;
                Atan = Sm1  - 2*mu*beta*Deviator ...
                       - 3*mu*(gamma-beta)*(sigdelas([1,2,4])*sigdelas([1,2,4])')/(sigeqelas)^2;
             end
             %% End of radial algo
%             if j==1 sigmaegal = sig( [3*j-2,3*j-1,3*j], i+1 )
%             sigmadevegal = sigd( [4*j-3,4*j-2,4*j-1,4*j], i+1 )
%             epsilonpegal = epsp( [3*j-2,3*j-1,3*j], i+1 ) end
             Ktane = S*Be'*Atan*Be;
             finte = S*Be'*sig( [3*j-2,3*j-1,3*j], i+1 );

             Ktan(map,map) = Ktan(map,map) + Ktane;
             fint(map) = fint(map) + finte;
          else
             error('Higher order elements unimplemented for elastoplasticity');
          end
       end
%       fint2 = fint(2)
       fi  = [ f(1:2*nnodes,i+1) ; zeros(nbloq,1) ] - fint - ...
             [ lambda(:,i+1) ; zeros(nbloq,1) ] ; % Don't re-imput the Dirichlet BC
       
       % Substract the internal forces to lambda ('cause it will be added later)
       lambda(:,i+1) = lambda(:,i+1) + Pro*fi(1:2*nnodes);
       
       % Project fint in order to suppress loading on the BCs
       %fi(1:2*nnodes) = fi(1:2*nnodes) - Pro*fi(1:2*nnodes);
       
       % End loop tests
       res = norm( fi(1:2*nnodes) );
       if res <= eps*resref
          stay = 0;
       end
       
       niter = niter+1;
       if niter > nmax
          stay = 0;
          warning([ 'Newton algo terminated without the desired convergence : ', num2str(res/resref) ])
       end
    end
    if talk == 1
       disp([ 'Time step ',num2str(i), ' finished with ',num2str(niter-1), ' elastoplastic iteration(s)' ]);
    end
    u(:,i+1)    = u(:,i) + Deltau;
    p(:,i+1)    = p(:,i) + Deltap;
    X(:,i+1)    = X(:,i) + DeltaX;
    
 end
 
end