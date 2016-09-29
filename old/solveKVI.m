function u = solveKVI(nodes,elements,E,nu,order,boundary,dirichlet,...
    tau,alpha,ntoelem,f,T,nt,varargin)

% This function solves the Kelvin-Voigt viscoélasticity problem
 
 global dotu
 global dotf

 mo = 0;
 if numel(varargin)>0
     mo = cell2mat(varargin(1));
 end
 
 nnodes = size(nodes,1);
 dt = T/nt;
 [K,C,nbloq,~,~] = Krig (nodes,elements,E,nu,order,boundary,...
     dirichlet,mo);
 
 u = zeros(2*nnodes+nbloq,nt);
 dotu = u;
 
 Kinter = K(1:2*nnodes,1:2*nnodes);
 finter = f(1:2*nnodes,:);
 
 % Estimation of the derivative of f :
 dotf = f-f;
 dotf(:,1) = (f(:,2)-f(:,1)) / dt;
 dotf(:,end) = (f(:,end)-f(:,end-1)) / dt;
 for i=2:nt
     dotf(:,i) = (f(:,i+1)-f(:,i-1))/(2*dt);
 end
 
 dotfinter = dotf(1:2*nnodes,:);
 
 % init : find dotu(1)
 RHS0 = (1+alpha)*finter(:,1) + tau*dotfinter(:,1) - alpha*Kinter*u(1:2*nnodes,1);
 LHS0 = Kinter*tau;
 % Add Lagrange stuff :
 RHS       = [RHS0;dotf(2*nnodes+1:end,1)];
 LHS       = [LHS0,C;C',zeros(nbloq)];
 dotu(:,1) = LHS\RHS;
 
 % Cranck-Nicholson scheme :
 for i=2:nt+1     
     RHS0 = (1+alpha)*finter(:,i) + tau*dotfinter(:,i) -...
         alpha*dt/2*Kinter*dotu(1:2*nnodes,i-1) - alpha*Kinter*u(1:2*nnodes,i-1);
     LHS0 = Kinter*(alpha*dt/2+tau);
     
     % Add Lagrange stuff :
     RHS = [RHS0;dotf(2*nnodes+1:end,i)];
     LHS = [LHS0,C;C',zeros(nbloq)];
     
     uin = LHS\RHS;
     dotu(1:2*nnodes,i) = uin(1:2*nnodes);
     u(:,i) = dt/2*(dotu(:,i)+dotu(:,i-1)) + u(:,i-1);
 end
 
end