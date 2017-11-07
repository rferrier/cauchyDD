% 26/10/2017
% Détection de fissure quelconque par écart à la réciprocité
% Régularisation par un algorithme génétique plus simple

tic
close all;
clear all;

warning('off','all'); % Because on matlab it's unbearable in picard stuff

addpath(genpath('./tools'))

% Parameters
jmax       = 0;     % Eigenvalues truncation number (if 0, automatic Picard choose)
ncrack     = 3;    % nb of cracks (odd : 1 crack, even : 2 cracks), 5 : 1% noise, 7 : 10% noise, 9 : corner crack, 11 : U crack
                    % 51  : refined postpro mesh
                    % 101 : direct & integrals refined, basic crack
                    % 103 : idem for the small crack
                    % 1001 : more test-functions
                    % 10001 : analysis on a mesh marking the crack
co         = [1,1,1,1]; % Coefficients for each RG test-case
niter      = 100;
popmin     = 40;  % Nb of survivors at each step (20,30)
multip     = 4;   % Population multiplicator at the mutation step (4,2)
crossing   = 1;   % Crossing multiplicator
pop        = popmin*(multip+1); % Population before the selection step
taumut     = .003;   % Mutation coefficient .003 = 1 elem only
Lregm1     = 1/20;%1/20; % Regularization length (1/20 mm-1)
mugrad     = 400;%100;%400;% Gradient penalization parameter (40 mm-1 or no unit (depends on gradnorm))
gradnorm   = 2;    % if gradnorm = 1 : Total Variation, 2 : squared, etc...
gradmode   = 0; % If gradmode == 1, use the normalized norm instead
gradmodemu = .0005; % Regularization parameter associated with gradmode = 1

%% In order to build the (Petrov-Galerkin) Left Hand Side
if ncrack == 10001
   [ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_ncu.msh' );
elseif ncrack < 50 || ncrack > 99
   [ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu.msh' );
else
   [ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu_r.msh' );
end
nnodesu = size(nodesu,1);

nelemu = size(elementsu,1); nn = size(elementsu,2); %nn=3

% Build the list of segment elements
nseg = nn*nelemu; segel = zeros(nseg,2); nbu = size(boundaryu,1);
for j=1:nelemu
   segel(nn*(j-1)+1,1) = elementsu(j,2);
   segel(nn*(j-1)+1,2) = elementsu(j,3);
   segel(nn*(j-1)+2,1) = elementsu(j,1);
   segel(nn*(j-1)+2,2) = elementsu(j,3);
   segel(nn*(j-1)+3,1) = elementsu(j,2);
   segel(nn*(j-1)+3,2) = elementsu(j,1);
end
j = 1;
while j <= nseg % Remove redundancy (I don't use unique because of inverted values)
   lis1 = find( segel(1:j-1,:) == segel(j,1) );
   lis1( find(lis1>j-1) ) = lis1( find(lis1>j-1) ) - j + 1;
   lis2 = find( segel(1:j-1,:) == segel(j,2) );
   lis2( find(lis2>j-1) ) = lis2( find(lis2>j-1) ) - j + 1;
   
   % also remove boundary elements
   lis3 = find( boundaryu(:,2:3) == segel(j,1) );
   lis3( find(lis3>nbu) ) = lis3( find(lis3>nbu) ) - nbu;
   lis4 = find( boundaryu(:,2:3) == segel(j,2) );
   lis4( find(lis4>nbu) ) = lis4( find(lis4>nbu) ) - nbu;
   
   if min(size( intersect(lis1,lis2) )) > 0 || min(size( intersect(lis3,lis4) )) > 0% || (size(lis3,1) > 0 || size(lis4,1) > 0)%
      segel(j,:) = [];
      j = j-1;
      nseg = nseg-1;
   end
   j = j+1;
end

% Build the shortcut list of ccordinates of segment elements
coorseg = [ nodesu(segel(:,1),1), nodesu(segel(:,1),2), ...
            nodesu(segel(:,2),1), nodesu(segel(:,2),2), ...
            (nodesu(segel(:,1),1)+nodesu(segel(:,2),1))/2, ...
            (nodesu(segel(:,1),2)+nodesu(segel(:,2),2))/2 ];

if ncrack == 1
   Anb = load('fields/matrix_20.mat'); %Anb = load('fields/matrix.mat');
   % The mesh is still needed
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
elseif ncrack == 2
   Anb = load('fields/matrix2_20.mat'); %Anb = load('fields/matrix2.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc2.msh' );
elseif ncrack == 3
   Anb = load('fields/matrix3_20.mat'); %Anb = load('fields/matrix3.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc3.msh' );
elseif ncrack == 5
   Anb = load('fields/matrix5_20.mat'); % Anb = load('fields/matrix5.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
elseif ncrack == 7
   Anb = load('fields/matrix7_20.mat'); %Anb = load('fields/matrix7.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
elseif ncrack == 9
   Anb = load('fields/matrix9_20.mat'); %Anb = load('fields/matrix9.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc9.msh' );
elseif ncrack == 11
   Anb = load('fields/matrix11_20.mat'); %Anb = load('fields/matrix11.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc11.msh' );
elseif ncrack == 51
   Anb = load('fields/matrix51.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
elseif ncrack == 101
   Anb = load('fields/matrix101.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r.msh' );
elseif ncrack == 103
   Anb = load('fields/matrix103.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_3.msh' );
elseif ncrack == 109
   Anb = load('fields/matrix109.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_9.msh' );
elseif ncrack == 111
   Anb = load('fields/matrix111.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_11.msh' );
elseif ncrack == 1001
   Anb = load('fields/matrix1001.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
elseif ncrack == 10001
   Anb = load('fields/matrix10001.mat');
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
end
Lhs = Anb.Lhs; Rhs1 = Anb.Rhs1; Rhs2 = Anb.Rhs2;
Rhs3 = Anb.Rhs3; Rhs4 = Anb.Rhs4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve the linear system and recover the crack
Rhs = Rhs1;
A = Lhs'*Lhs; b = Lhs'*Rhs; szebai = size(b,1);
b1 = Lhs'*Rhs1; b2 = Lhs'*Rhs2; b3 = Lhs'*Rhs3; b4 = Lhs'*Rhs4;

%% Build the regularization matrix
% Initialize the vectors
nogap1  = ones(nseg,1); nogap2  = ones(nseg,1);
nogap3  = ones(nseg,1); nogap4  = ones(nseg,1);

oldauthorized = cell(pop,1); oldauthorized2 = cell(pop,1);
nodescrack    = cell(pop,1); nodescrack2    = cell(pop,1);

nogap1 = zeros(nseg,pop); nogap2 = zeros(nseg,pop);
nogap3 = zeros(nseg,pop); nogap4 = zeros(nseg,pop);

res1 = zeros(pop,1); res2 = zeros(pop,1); % Cost function stuff
res3 = zeros(pop,1); res4 = zeros(pop,1);
res0 = zeros(pop,1); lres = zeros(pop,1);
ressum = zeros(pop,1); gradsum = zeros(pop,1);
gradres1 = zeros(pop,1); gradres2 = zeros(pop,1);
gradres3 = zeros(pop,1); gradres4 = zeros(pop,1);

% Compute lengths and normals
lengvect1 = zeros(nseg,1);
normalct1 = zeros(nseg,2);
for k=1:nseg
   no1 = segel(k,1); no2 = segel(k,2);
   lengvect1(k) = sqrt( ( nodesu(no1,1)-nodesu(no2,1) )^2 +...
                       ( nodesu(no1,2)-nodesu(no2,2) )^2 );
   normalct1(k,:) = [ nodesu(no1,2)-nodesu(no2,2) ,...
                    -(nodesu(no1,1)-nodesu(no2,1)) ] / lengvect1(k) ;
end

for j=1:pop % Randomly initialize the population
   thisseg = 1+floor(nseg*rand);
   oldauthorized{j} = thisseg;
   nodescrack{j}    = [segel(thisseg,1);segel(thisseg,2)];
end

tic
for i=1:niter
   %% Selection stage : compute the criterion
   nogap1 = zeros(nseg,pop); nogap2 = zeros(nseg,pop);
   nogap3 = zeros(nseg,pop); nogap4 = zeros(nseg,pop);
   
   for j=1:pop
      oldauthorized0 = unique(oldauthorized{j}); % Just for safety
      authorized = [ 2*oldauthorized0-1 ; 2*oldauthorized0 ];% authorized = authorized(:);
%      nodescrack0 = unique(nodescrack{j});
      nodescrack0 = unique( [segel(oldauthorized0,1);segel(oldauthorized0,2)] );

      %% compute the criterion for each of them
      [Q,Theta] = eig( A(authorized,authorized) );
      thetas = diag(Theta);
      [thetas,Ind] = sort( thetas,'descend' );
      Q = Q(:,Ind);
      Thetas = diag(thetas); 
      Theta = Theta(Ind,Ind);

      if size(authorized,1)==0 || size(authorized,2)==0 % Exception : no element at all
         res1(j) = norm( b1 ); res2(j) = norm( b2 );
         res3(j) = norm( b3 ); res4(j) = norm( b4 );
         lres(j) = 0;
         res0(j) = (res1(j)*co(1) + res2(j)*co(2) + res3(j)*co(3) + res4(j)*co(4)) / ...
                   (norm(b1)*co(1) + norm(b2)*co(2) + norm(b3)*co(3) + norm(b4)*co(4));
      else
         % Compute Picard stopping indices
         if size(authorized,1)==2
            ind1 = 2; ind2 = 2; ind3 = 2; ind4 = 2;
         else % If there are more than 2 dofs, use Picard
      
            % Plot the Picard stuff
            imax = min( find(thetas/thetas(1)<1e-16) );
            if size(imax,1) == 0
               imax = size(thetas,1);
            end

            tplo = thetas(1:imax); bplo1 = Q'*b1(authorized);
            rplo1 = (bplo1)./thetas; bplo1 = bplo1(1:imax);
            rplo1 = rplo1(1:imax); bplo2 = Q'*b2(authorized);
            rplo2 = (bplo2)./thetas; bplo2 = bplo2(1:imax);
            rplo2 = rplo2(1:imax); bplo3 = Q'*b3(authorized);
            rplo3 = (bplo3)./thetas; bplo3 = bplo3(1:imax);
            rplo3 = rplo3(1:imax); bplo4 = Q'*b4(authorized);
            rplo4 = (bplo4)./thetas; bplo4 = bplo4(1:imax);
            rplo4 = rplo4(1:imax);

            % Remove Zeros in rploi (why on hell are there zeros in the first place ?)
            me1 = mean(abs(rplo1))/1e5; arplo1 = max(me1,abs(rplo1));
            me2 = mean(abs(rplo2))/1e5; arplo2 = max(me2,abs(rplo2));
            me3 = mean(abs(rplo3))/1e5; arplo3 = max(me3,abs(rplo3));
            me4 = mean(abs(rplo4))/1e5; arplo4 = max(me4,abs(rplo4));

            ind1 = findPicard2 (log10(arplo1), 7, 1, 3);
            ind2 = findPicard2 (log10(arplo2), 7, 1, 3);
            ind3 = findPicard2 (log10(arplo3), 7, 1, 3);
            ind4 = findPicard2 (log10(arplo4), 7, 1, 3);
         end
      
         % Filter eigenvalues
         if jmax == 0
           jmax0 = size(Thetas,1);
            jmax1 = ind1; jmax2 = ind2; jmax3 = ind3; jmax4 = ind4;
%            jmax1 = jmax0; jmax2 = jmax0; jmax3 = jmax0; jmax4 = jmax0;
         else
            jmax0 = min( size(Thetas,1) , jmax );
            jmax1 = jmax; jmax2 = jmax; jmax3 = jmax; jmax4 = jmax;
         end
      
         %ThetaT  = Thetas( 1:jmax0 , 1:jmax0 );
         ThetaT1 = Thetas( 1:jmax1 , 1:jmax1 );
         ThetaT2 = Thetas( 1:jmax2 , 1:jmax2 );
         ThetaT3 = Thetas( 1:jmax3 , 1:jmax3 );
         ThetaT4 = Thetas( 1:jmax4 , 1:jmax4 );
      
         bT1 = Q'*b1(authorized); bT1 = bT1(1:jmax1);
         bT2 = Q'*b2(authorized); bT2 = bT2(1:jmax2);
         bT3 = Q'*b3(authorized); bT3 = bT3(1:jmax3);
         bT4 = Q'*b4(authorized); bT4 = bT4(1:jmax4);
      
         Solu1 = zeros(szebai,1); Solu2 = zeros(szebai,1);
         Solu3 = zeros(szebai,1); Solu4 = zeros(szebai,1);

         Solu1(authorized) = Q(:,1:jmax1) * (ThetaT1\bT1);
         Solu2(authorized) = Q(:,1:jmax2) * (ThetaT2\bT2);
         Solu3(authorized) = Q(:,1:jmax3) * (ThetaT3\bT3);
         Solu4(authorized) = Q(:,1:jmax4) * (ThetaT4\bT4);

         res1(j) = norm( A*Solu1 - b1 );% / norm( b1 );
         res2(j) = norm( A*Solu2 - b2 );% / norm( b2 );
         res3(j) = norm( A*Solu3 - b3 );% / norm( b3 );
         res4(j) = norm( A*Solu4 - b4 );% / norm( b4 );

         % Compute length of the crack
         lres(j) = 0;
         for k=1:size(oldauthorized0)
            lres(j) = lres(j) + lengvect1(oldauthorized0(k));
         end

         % Compute gradient cost function (variance)
         gradres1(j) = 0; gradres2(j) = 0; gradres3(j) = 0; gradres4(j) = 0;
         for k=1:size(nodescrack0,1)
            no1 = nodescrack0(k);
            segs0 = rem( find(segel(oldauthorized0,:)==no1)-1 , size(oldauthorized0,1) )+1;
            segs = oldauthorized0(segs0);
            if size(segs,1) == 1
               gradloc1 = ([Solu1(2*segs-1);Solu1(2*segs)]'*normalct1(segs,:)')/lengvect1(segs);
               gradloc2 = ([Solu2(2*segs-1);Solu2(2*segs)]'*normalct1(segs,:)')/lengvect1(segs);
               gradloc3 = ([Solu3(2*segs-1);Solu3(2*segs)]'*normalct1(segs,:)')/lengvect1(segs);
               gradloc4 = ([Solu4(2*segs-1);Solu4(2*segs)]'*normalct1(segs,:)')/lengvect1(segs);
            else
               NormVal1 = diag([Solu1(2*segs-1),Solu1(2*segs)]*normalct1(segs,:)');
               NormVal2 = diag([Solu2(2*segs-1),Solu2(2*segs)]*normalct1(segs,:)');
               NormVal3 = diag([Solu3(2*segs-1),Solu3(2*segs)]*normalct1(segs,:)');
               NormVal4 = diag([Solu4(2*segs-1),Solu4(2*segs)]*normalct1(segs,:)');
               gradloc1 = var(NormVal1)/mean(lengvect1(segs))*sqrt(size(segs,1));
               gradloc2 = var(NormVal2)/mean(lengvect1(segs))*sqrt(size(segs,1));
               gradloc3 = var(NormVal3)/mean(lengvect1(segs))*sqrt(size(segs,1));
               gradloc4 = var(NormVal4)/mean(lengvect1(segs))*sqrt(size(segs,1));
            end

            gradres1(j) = gradres1(j) + abs(gradloc1)^gradnorm*mean(lengvect1(segs));
            gradres2(j) = gradres2(j) + abs(gradloc2)^gradnorm*mean(lengvect1(segs));
            gradres3(j) = gradres3(j) + abs(gradloc3)^gradnorm*mean(lengvect1(segs));
            gradres4(j) = gradres4(j) + abs(gradloc4)^gradnorm*mean(lengvect1(segs));
         end

         gradsum(j) = (gradres1(j)*co(1) + gradres2(j)*co(2) + gradres3(j)*co(3) + gradres4(j)*co(4)) / (sum(co));
         ressum(j) = (res1(j)*co(1) + res2(j)*co(2) + res3(j)*co(3) + res4(j)*co(4)) / ...
                     (norm(b1)*co(1) + norm(b2)*co(2) + norm(b3)*co(3) + norm(b4)*co(4));

         res0(j) = mugrad*gradsum(j) + ressum(j) + lres(j)*Lregm1;

        % Re-loop over segments to build gap's norm
         for kk = 1:size(oldauthorized0,1)
            k = oldauthorized0(kk);
            no1 = segel(k,1); no2 = segel(k,2);
            nogap1(k,j) = norm( Solu1( [2*k-1,2*k] ) );
            nogap2(k,j) = norm( Solu2( [2*k-1,2*k] ) );
            nogap3(k,j) = norm( Solu3( [2*k-1,2*k] ) );
            nogap4(k,j) = norm( Solu4( [2*k-1,2*k] ) );
         end
      end
   end

   %% Selection stage : kill the less efficient (rem : less efficient will be overwritten later)
   % Reorder the coefs
   [res0,Ind] = sort( res0 );
   res1 = res1(Ind); res2 = res2(Ind); res3 = res3(Ind); res4 = res4(Ind);
   gradres1 = gradres1(Ind); gradres2 = gradres2(Ind);
   gradres3 = gradres3(Ind); gradres4 = gradres4(Ind);
   ressum = ressum(Ind); lres = lres(Ind); gradsum = gradsum(Ind);
   nogap1 = nogap1(:,Ind); nogap2 = nogap2(:,Ind);
   nogap3 = nogap3(:,Ind); nogap4 = nogap4(:,Ind);

   for j=1:size(Ind) % It doesn't work in a raw...
     oldauthorized2{j} = oldauthorized{Ind(j)};
     nodescrack2{j}    = nodescrack{Ind(j)};
   end
   oldauthorized = oldauthorized2; nodescrack = nodescrack2;
   
   %% Mutation Stage
   rank = popmin+1; % Rank of the current individual to add

   for j=1:popmin
      for k=1:multip
         oldauthorized0 = unique(oldauthorized{j});
         tochange = rand(floor(taumut*nseg),1);
         tochange = floor(tochange*nseg)+1;

         % Change the prescripted bits
         tochange2      = setdiff(tochange,oldauthorized0);
         oldauthorized0 = setdiff(oldauthorized0,tochange);
         oldauthorized0 = union(oldauthorized0,tochange2);
         if isrow(oldauthorized0)
            oldauthorized0 = oldauthorized0';
         end
         oldauthorized{rank} = unique(oldauthorized0);
         rank = rank+1;
      end
   end

   % Remove double and add monoms in their place (unique doesn't work)
%   zzz = 0;
   for j=2:pop
      oldauthorizedj = oldauthorized{j};
      for k=1:j-1
         oldauthorizedk = oldauthorized{k};
         if size(oldauthorizedj,1) == size(oldauthorizedk,1) && size(oldauthorizedj,2) == size(oldauthorizedk,2)
            if sum( oldauthorizedj==oldauthorizedk ) == size(oldauthorizedj,1)
               % Suppress this double (put a fresh new one on it)
               thisseg = 1+floor(nseg*rand);
               oldauthorized{j} = thisseg;
               nodescrack{j}    = [segel(thisseg,1);segel(thisseg,2)];
%               zzz++;
               break;
            end
         end
      end
   end
end
disp([ 'Genetic algo terminated ', num2str(toc) ]);
% Choose the best one
%[minresid, chosen] = min(res0);
chosen = 1;
oldauthorized0 = oldauthorized{chosen}; 
nogap10 = nogap1(:,chosen); nogap20 = nogap2(:,chosen);
nogap30 = nogap3(:,chosen); nogap40 = nogap4(:,chosen);

% Segments visu

figure;
hold on;
plot( [nodes(5,1),nodes(6,1)], [nodes(5,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
if ncrack == 2
   plot( [nodes(7,1),nodes(8,1)], [nodes(7,2),nodes(8,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
elseif ncrack == 9 || ncrack == 109 || ncrack == 11 || ncrack == 111
   plot( [nodes(7,1),nodes(6,1)], [nodes(7,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
end
if ncrack == 11 || ncrack == 111
   plot( [nodes(7,1),nodes(8,1)], [nodes(7,2),nodes(8,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
end
nogapp1 = co(1)*nogap10 + co(2)*nogap20 + co(3)*nogap30 + co(4)*nogap40;
maxn1 = max(nogapp1);
%      for i=1:nseg
for j=1:size(oldauthorized0,1)
   i = oldauthorized0(j);
   no1 = segel(i,1); no2 = segel(i,2);
   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
   
   x = nogapp1(i)/maxn1;
   rgb = rgbmap(x);
   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
end
%axis equal;
axis([0 1 0 1]);
colormap('default')
h = colorbar();
ytick = get (h, 'ytick');
set (h, 'yticklabel', sprintf ( '%g|', maxn1*ytick+min(nogap10) ));

%total = res0(1)
%gradient = mugrad*gradsum(1)
%no_reg = ressum(1)
%length = lres(1)*Lregm1

%total = res0(1)
%gradient = gradmodemu*gradsum(1)
%no_reg = ressum(1)

%figure;
%hold on;
%plot(res1);
%plot(res2,'Color','black');
%
%figure;
%hold on;
%plot(ind1);
%plot(ind2,'Color','black');
%plot(ind3,'Color','red');
%plot(ind4,'Color','green');
%legend('Picard index 1','Picard index 2','Picard index 3','Picard index 4');