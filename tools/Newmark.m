function u = Newmark (M, C, K, f, u0, v0, a0, dt, beta, gamma)
% This function solves the dynamic system with a Newmark approach
% Ma + Cv + Ku = f

% input : f     : mulit-row matrix of the loading
%         M     : mass matrix
%         C     : Damping matrix
%         K     : stiffness matrix (with its Lagrange multipliers)
%         u0    : displacement at t = 0
%         v0    : velocity at t = 0
%         a0    : accelerations at t = 0
%         dt    : time discretization parameter (constant)
%         beta  : parameter of the method
%         gamma : parameter of the method
 
% output : u    : multi-row matrix of the solution

 ndof  = size(f,1);
 ntime = size(f,2);

 u = zeros( ndof, ntime );
 up = u0; vp = v0; ap = a0;
 
 fatK = K + 1/(beta*dt^2)*M + gamma/(beta*dt)*C;
 
 for i=1:ntime
    u(:,i) = fatK \ ( f(:,i) + ...
               C*( gamma/(beta*dt)*up + (gamma/beta-1)*vp + dt/2*(gamma/beta-1)*ap ) + ...
               M*( 1/(beta*dt^2)*up + 1/(beta*dt)*vp + (1/(2*beta)-1)*ap ) );
               
    ac = 1/(2*beta*dt^2/2) * ( u(:,i) - up - dt*vp - dt^2/2*(1-2*beta)*ap );
    vc = vp + dt*( (1-gamma)*ap + gamma*ac );
    
    up = u(:,i); vp = vc; ap = ac;
 end

end
