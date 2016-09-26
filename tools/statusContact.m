function [u] = statusContact( nodes, elem, mat, order, boundary,...
    dirichlet, eta, nmax, ntoelem, fneumann, varargin )
% This function solves the contact problem with the status method
% dirichlet has a new format :
%    first index  : index of the boundary
%    second index : 1=x, 2=y
%    third        : value
%    fourth       : 0 : equal, 1 : superior, -1 : inferior

 mo = 0;
 if numel(varargin)>0
     mo = cell2mat(varargin(1));
 end

nnodes = size(nodes,1);

%% First find witch node must get Dirichlet
dofD     = zeros(2*nnodes,4);
cont2dof = [];% Contains the indices of the Dirichlet boundaries
n = 0;
for i = 1:size(boundary,1)
   for j=1:size(dirichlet,1)
      if boundary(i,1) == dirichlet(j,1)
         no1 = boundary(i,2);
         no2 = boundary(i,3);
         if dirichlet(j,2) == 1 % x
            dofD(2*no1-1,1) = 1;
            dofD(2*no1-1,2) = dirichlet(j,3); % Value
            dofD(2*no1-1,3) = dirichlet(j,4); % > < =
            if size(find(cont2dof == 2*no1-1),2) == 0
               n = n+1;
               cont2dof(n) = 2*no1-1;
               dofD(2*no1-1,4) = n;
            end
            
            dofD(2*no2-1,1) = 1;
            dofD(2*no2-1,2) = dirichlet(j,3); % Value
            dofD(2*no2-1,3) = dirichlet(j,4); % > < =
            if size(find(cont2dof == 2*no2-1),2) == 0
               n = n+1;
               cont2dof(n) = 2*no2-1;
               dofD(2*no2-1,4) = n;
            end
            
         elseif dirichlet(j,2) == 2 % y
            dofD(2*no1,1) = 1;
            dofD(2*no1,2) = dirichlet(j,3); % Value
            dofD(2*no1,3) = dirichlet(j,4); % > < =
            if size(find(cont2dof == 2*no1),2) == 0
               n = n+1;
               cont2dof(n) = 2*no1;
               dofD(2*no1,4) = n;
            end
            
            dofD(2*no2,1) = 1;
            dofD(2*no2,2) = dirichlet(j,3); % Value
            dofD(2*no2,3) = dirichlet(j,4); % > < =
            if size(find(cont2dof == 2*no2),2) == 0
               n = n+1;
               cont2dof(n) = 2*no2;
               dofD(2*no2,4) = n;
            end
            
         end
      end
   end
end

% Store the nodes with Statuts to do : don't change this matrix.
dofC = dofD;

%% In case friction is activated : block every potentially blocked dof.
% /!\ Problem if there is a Dirichlet on the other side
if eta ~= 0
   for j=1:2*nnodes
      if dofD(j,3) ~= 0
         m = j + 2*mod(j,2) - 1; % the other component
         cont2dof( end+1 ) = m;
         dofD(m,:) = [1,0,0,n+1];  % Impose ths one
         n = n+1;
      end
   end
end

%% Loop over Status iterations
nochange = 0; % Just in case nmax<1
fcoulomb = zeros(2*nnodes,1);
for i=1:nmax
   % Build C
   C = zeros(2*nnodes, n);
   for j=1:n
      C(cont2dof(j), j) = 1;
   end
   
   % Build Dirichlet RHS
   fdir = zeros(2*nnodes + n,1);
   for j=1:n
      fdir(2*nnodes+j) = dofD(cont2dof(j), 2);
   end
   
   % Build K
   [K,~,~,~,~] = Krig2 (nodes, elem, mat, order, boundary, [], mo);
   
   % Assembly stuff
   K = [K,C ; C',zeros(n)];
   
   % Solve the problem
   u = K \ ( fdir + [fneumann;zeros(n,1)] + [fcoulomb;zeros(n,1)]);

   % Extract RHS :
   frea = zeros(2*nnodes,1);
   for j=1:n
      frea(cont2dof(j)) = -u(2*nnodes+j);
   end
   
   nochange  = 1;
   fcoulombo = fcoulomb;
   fcoulomb  = zeros(2*nnodes,1);
   % Test the status
   for j=1:2*nnodes
      if dofC(j,3) < 0 % <
         if frea(j) > 0 && dofD(j,1) == 1
            cont2dof( dofD(j,4) ) = [];
            % Renumber dofD
            for k=1:2*nnodes
               if dofD(k,4) > dofD(j,4)
                  dofD(k,4) = dofD(k,4)-1;
               end
            end
            dofD(j,:) = 0;  % Don't impose ths one
            n = n-1;
            nochange = 0;
         elseif u(j) > dofC(j,2) && dofD(j,1) == 0
            cont2dof( end+1 ) = j;
            dofD(j,:) = dofC(j,:);  % Impose ths one
            n = n+1;
            dofD(j,4) = n;
            nochange = 0;
         end
      elseif dofC(j,3) > 0
         if frea(j) < 0 && dofD(j,1) == 1
            cont2dof( dofD(j,4) ) = [];
            % Renumber dofD
            for k=1:2*nnodes
               if dofD(k,4) > dofD(j,4)
                  dofD(k,4) = dofD(k,4)-1;
               end
            end
            dofD(j,:) = 0;  % Don't impose ths one
            n = n-1;
            nochange = 0;
         elseif u(j) < dofC(j,2) && dofD(j,1) == 0
            cont2dof( end+1 ) = j;
            dofD(j,:) = dofC(j,:);  % Impose ths one
            n = n+1;
            dofD(j,4) = n;
            nochange = 0;
         end
      end
      
      % Coulomb stuff
      if dofD(j,3) ~= 0
         % Identifie the other component
         m = j + 2*mod(j,2) - 1;
         if abs(frea(m)) > abs(eta*frea(j)) && dofD(m,1) == 1  % /!\ Problem if k is also Dirichletized
            cont2dof( dofD(m,4) ) = [];
            % Renumber dofD
            for k=1:2*nnodes
               if dofD(k,4) > dofD(m,4)
                  dofD(k,4) = dofD(k,4)-1;
               end
            end
            dofD(m,:) = 0;  % Don't impose ths one
            n = n-1;
            nochange = 0;
            fcoulomb(m) = frea(m)*eta*abs(frea(j)/frea(m));
         elseif u(m)*fcoulombo(m) > 0 && dofD(m,1) == 0 % Work must be < 0
            cont2dof( end+1 ) = m;
            dofD(m,:) = [1,0,0,n+1];  % Impose ths one
            n = n+1;
            nochange = 0;
         end
      end
   end
   disp([ 'Contact Iteration ',num2str(i), ' finished.' ]);

   % Test end of algo
   if nochange
      break;
   end
   
end

if nochange ~= 1
   warning('Status algorithm ended before convergence')
end

end