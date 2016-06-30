function u = solveDam( nodes, elements, E, nu, order, boundary,...
    dirichlet, neumann, damage, error, nmax )
% This function solves the non-linear problem corresponding to a damaged
% composite material : young modulus is stronger in compression than in
% traction.

 % The right hand side and first step :
 [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
 f = loading(nbloq,nodes,boundary,neumann);
 u = K\f;

 % Estimate the residual : 
[K,C,nbloq] = KrigDam (nodes,elements,E,nu,order,boundary,...
	dirichlet,u,damage);
 res = f - K*u;

 iter = 0;
 while norm(res) > error
     if iter > nmax
         warning('Non-linear resolution ended before convergence !')
         break
     end
     u = u + K\res;
     
     % Estimate the residual : 
     [K,C,nbloq] = KrigDam (nodes,elements,E,nu,order,boundary,...
	     dirichlet,u,damage);
     res = f - K*u;
     
     iter = iter+1;
 end
iter
end
