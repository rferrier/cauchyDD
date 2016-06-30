function u = solveHyp( nodes, elements, E, nu, order, boundary,...
    dirichlet, alpha, error, nmax, ntoelem, fdir, fneumann )
% This function solves the non-linear problem corresponding to an
% hyperelastic material with compressibility nonlinearity
% /!\ Plane deformations

 nnodes = size(nodes,1);

 % The right hand side and first step :
 [K,C,~] = Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
 u = K\(fneumann+fdir);

 % Estimate the residual :
 lagr = C*u(2*nnodes+1:end);
 rhs = stressHyp( u,E,nu,nodes,elements,order,2,ntoelem,alpha,1 ) + lagr;
 res = fneumann(1:2*nnodes) - rhs;
 
%  hold on
%  plot(res)
%  plot(lagr,'Color','red')
%  plot(f(1:2*nnodes),'Color','black')
 
 err = zeros(nmax+1,1);
 err(1) = norm(res)/norm( ( fneumann(1:2*nnodes) + lagr ) );
 
 iter = 0;
 while err(iter+1) > error
     if iter >= nmax
         warning('Non-linear resolution ended before convergence !')
         break
     end
     
     [K,C,nbloq] = KrigHyp (nodes,elements,E,nu,order,boundary,...
         dirichlet,alpha, u, 1);
     u = u + K \ [res ; zeros(nbloq,1) ];
     
     % Estimate the residual :
     lagr = C*u(2*nnodes+1:end);
     rhs = stressHyp( u,E,nu,nodes,elements,order,2,ntoelem,alpha,1 ) + lagr;
     res = fneumann(1:2*nnodes) - rhs;
     err(iter+2) = norm(res)/norm( ( fneumann(1:2*nnodes) + lagr ) );
     iter = iter+1;
     
 end
 %iter
 %plot(log(err))
end
