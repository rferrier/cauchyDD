function [u] = statusIneq( Ki, C, b, f, nmax )
% This function solves the inequality BC problem with the status method
% No friction for now.
% Ki already contains hard Dirichlet BC
% C and b are such that C'u <= b

 nnodes = size(Ki,1);  % Year, that's not exactly nnodes, but...
 ncont  = size(C,2);

 n       = ncont; % init of n
 statusL = ones(ncont,1); % Stores the status

 %% Loop over Status iterations
 nochange = 0; % Just in case nmax<1
 for i=1:nmax
   
    % Build Cc and Dirichlet RHS
    Cc = C; bc = b;
    torm = []; % Lines to remove in Cc and bc
    for j=1:ncont
       if statusL(j) == 0
          torm(end+1,1) = j;
       end
    end
    Cc(:,torm) = []; bc(torm) = [];

    % Assembly stuff
    K = [Ki,Cc ; Cc',zeros(n)];
   
    % Solve the problem
    u = K \ ( [zeros(nnodes,1);bc] + [f;zeros(n,1)] );

    % Extract diverse things :
    uint = u( 1:nnodes );
    frea = u( nnodes+1:end );

    nochange = 1;
    
    % Test the status
    freac = C'*Cc*frea;  % freac > 0
    impdof = C'*uint;       %impdof <= b
   
    for j=1:ncont
       if freac(j) < 0 && statusL(j,1) == 1
          statusL(j,1) = 0;  % Don't impose ths one
          n = n-1;
          nochange = 0;
       elseif impdof(j) > b(j) && statusL(j,1) == 0
          statusL(j,1) = 1;   % Impose ths one
          n = n+1;
          nochange = 0;
       end
   end
      
   %disp([ 'Status Iteration ',num2str(i), ' finished.' ]);
   
   % Test end of algo
   if nochange
      break;
   end
   
 end

 if nochange ~= 1
    warning('Status algorithm ended before convergence')
 end
end