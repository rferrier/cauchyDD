function [u, lagr] = solveDN(  )
 % Function that solves a particuliar DN problem
 f1in = loading(nbloq1,nodes,boundary,neumann1);
 fdir = dirichletRhs(u2, 3, C1, boundary);
 f1 = f1in + fdir;

 uin1 = K1\f1;
 u = uin1(1:2*nnodes,1);
 lagr = uin1(2*nnodes+1:end,1);

end

