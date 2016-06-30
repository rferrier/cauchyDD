function u = solveKVI(nodes,elements,E,nu,order,boundary,dirichlet,...
    tau,alpha,ntoelem,f,T,nt,varargin)

% This function solves the Kelvin-Voigt viscoélasticity vith Instantaneus
% elasticity problem
 
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
 RHS0 = (1+alpha)/tau*finter(:,1) + dotfinter(:,1) - alpha/tau*Kinter*u(1:2*nnodes,1);
 LHS = K;
 % Add Lagrange stuff :
 RHS = [RHS0;dotf(2*nnodes+1:end,1)];
 %RHS = (1+alpha)/tau*f(:,1) + dotf(:,1) - alpha/tau*K*u(:,1);
%  LHS       = [LHS0,tau*C;tau*C',zeros(nbloq)];
 dotu(:,1) = LHS\RHS;
 
 % Cranck-Nicholson scheme :
 for i=2:nt+1     
     RHS0 = (1+alpha)/(alpha*dt/2+tau)*finter(:,i) + tau/(alpha*dt/2+tau)*dotfinter(:,i) -...
         alpha/(alpha*dt/2+tau)*dt/2*Kinter*dotu(1:2*nnodes,i-1) - alpha/(alpha*dt/2+tau)*Kinter*u(1:2*nnodes,i-1);
     LHS = K;

%      % Add Lagrange stuff :
     RHS = [ RHS0 ; dotf(2*nnodes+1:end,i) ];% - alpha/(alpha*dt/2+tau)*dt/2*dotu(2*nnodes+1:end,i-1) ] ;% + (1+alpha)/(alpha*dt/2+tau)*f(2*nnodes+1:end,i) ];
     %RHS = (1+alpha)/(alpha*dt/2+tau)*f(:,i) + tau/(alpha*dt/2+tau)*dotf(:,i) -...
      %   alpha/(alpha*dt/2+tau)*dt/2*K*dotu(:,i-1) - alpha/(alpha*dt/2+tau)*K*u(:,i-1);
%      LHS = [LHS0,(tau+alpha*dt/2)*C;(tau+alpha*dt/2)*C',zeros(nbloq)];
     
     dotu(:,i) = LHS\RHS;
     u(:,i) = dt/2*(dotu(:,i)+dotu(:,i-1)) + u(:,i-1);
 end
 
end