function [ Y, Theta, flag ] = gradEig ( A, L, n, ker, varargin )
% This function computes the n first (generalized) eigenvalues of a matrix
% via a congugate gradient and Ritz analysis.
% Works well for n<<sA

% @param : A   : matrix to find the eigenelemnts
%          L   : Preconditionner
%          n   : size of the Ritz basis
%          ker : wether to test the kernel of L (time costly)

 sA = size(A,1);
 if size(A,2) ~= sA
    error("The matrix is rectangular")
 end

 flag = 0; % Means ended normally
 
 b  = rand(sA,1);   % Random Rhs
 if numel(varargin)>0
    b = cell2mat(varargin(1));
 end
 
 invLgiven = 0;
 if numel(varargin)>1
    invLgiven = cell2mat(varargin(2));
 end
 
 z = zeros(sA,n+1);
 r = zeros(sA,n+1);
 d = zeros(sA,n+1);
 
 pinvert = 0;
 if ker == 1 && invLgiven == 0
    %kL = fatNull(L); %PL = eye(sA) - kL*((kL'*kL)\kL');
    if rcond(L) < 1e-15 %size(kL,2) > 0
       pinvert = 1;
       pinvL = pinv(L);
    end
 end

 alpha = zeros(n+1,1); beta = zeros(n+1,1);
 
 x      = zeros(sA,1);  % Initialize solution
 r(:,1) = b-A*x;        % Residual
 
 if invLgiven == 1      % Actually, the invert of L was given
    z(:,1) = L*r(:,1);
 elseif pinvert == 0
    z(:,1) = L\r(:,1);     % Preconditionned residual
 else
    z(:,1) = pinvL*r(:,1);
 end
 
 d(:,1) = z(:,1);       % Search direction

 V = zeros(sA,n); H = zeros(n);
 eta = 0;
 
 for i=1:n
    Ad       = A*d(:,i);
    den      = d(:,i)'*Ad;
    
    if den < 0 % Floating point problem
       flag = 1; % Cannot make this iteration
       break;
    end
    
    sden     = sqrt(den);
    d(:,i)   = d(:,i)/sden;
    Ad       = Ad/sden;
    num      = r(:,i)'*d(:,i);
    x        = x + d(:,i)*num;
    r(:,i+1) = r(:,i) - Ad*num; 

    if invLgiven == 1      % Actually, the invert of L was given
       z(:,i+1) = L*r(:,i+1);
    elseif pinvert == 0
       z(:,i+1) = L\r(:,i+1);     % Preconditionned residual
    else
       z(:,i+1) = pinvL*r(:,i+1);
    end
    
    alpha(i) = num/sqrt(den);
    beta(i)  = - z(:,i+1)'*Ad/sden;
    
    % Reorthogonalization of residual
    for j=1:i
        betac = z(:,i+1)'*r(:,j) / (z(:,j)'*r(:,j));
        z(:,i+1) = z(:,i+1) - z(:,j) * betac;
        r(:,i+1) = r(:,i+1) - r(:,j) * betac;
    end
    
    %% Orthogonalization
    d(:,i+1) = z(:,i+1);
    for j=i:i % Mm...
        betaij = z(:,i+1)'*Ad;
        d(:,i+1) = d(:,i+1) - d(:,j) * betaij;
    end
    
    V(:,i) = (-1)^(i-1)*z(:,i)/(sqrt(r(:,i)'*z(:,i)));

    etap  = eta;
    delta = 1/alpha(i);
    if i > 1
       delta  = delta + beta(i-1)/alpha(i-1);
    end
    eta = sqrt(beta(i))/alpha(i);
    
    if i > 1
       H(i,[i-1,i]) = [etap, delta];
       H(i-1,i)     = etap;
    else
       H(i,i)       = delta;
    end
    
 end

 if flag == 1
    V(:,i:end) = [];
    H(:,i:end) = [];
    H(i:end,:) = [];
 end
 
  % Compute eigenelems of the Hessenberg
 [Q,Theta] = eig(H);
 theta = diag(Theta);
 [theta,Ind] = sort(theta,'descend');
 Q = Q(:,Ind);
 Theta = Theta(Ind,Ind);
 Y = V*Q;
 
% norm(A*x-b)
% norm(A*x)
% norm(r(:,end))
 
end