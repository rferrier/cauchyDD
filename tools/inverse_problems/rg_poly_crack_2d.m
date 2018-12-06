function [normal,Cte1,Cte2,ug,x] = rg_poly_crack_2d( nodes, extnorm, boundary, order, mat, u, f, ordp, v, varargin )
% This function identifies a segment crack with RG polynomial galerkin from 2 test-cases

% input : nodes    : nodes of the underlying mesh
%         extnorm  : exterior normal per boundary element 
%         boundary : boundaries of the underlying mesh
%         order    : order of the underlying mesh
%         mat      : used material
%         u        : Dirichlet data (x2)
%         f        : Neumann data (x2) given at each boundary element (not a nodal force, but lineic force, constant per element)
%         ordp     : nb of polynomial basis functions
%         v        : should the algorithm talk ?
%         varargin : if 1 : use Fourier instead of polynoms
%                    1 or 2 : test-case for the constant

% output : normal  : normal to the crack
%          Cte1    : cte giving the line from test-case 1
%          Cte2    : cte giving the line from test-case 2
%          ug      : gap of u on the line
%          x       : abscissa on the line (to plot ug)

 if size(u,2)~=2
    warning('u should be of size n*2');
 end

 if size(f,2)~=2
    warning('f should be of size n*2');
 end

 fourier = 0;
 if numel(varargin)>0
    fourier = cell2mat(varargin(1));
 end

 cst = 2;
 if numel(varargin)>1
    cst = cell2mat(varargin(2));
 end

 nboun = size(boundary,1);

 E = mat(2);
 nu = mat(3);

 nnodes = size(nodes,1);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% First task : find the normal

 if order==1
    Ng = 1;
 elseif order==2
    Ng  = 2; 
 end

 R11 = 0; R21 = 0; R31 = 0; R12 = 0; R22 = 0; R32 = 0;
 for i=1:nboun
    bonod = boundary(i,:); exno = extnorm(i,:)';
 
    no1 = bonod(2); no2 = bonod(3);
    x1 = nodes(no1,1); y1 = nodes(no1,2);
    x2 = nodes(no2,1); y2 = nodes(no2,2);
    len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
    
 %   % Compute sigma.n
    nomat = [ exno(1), 0, exno(2) ; ...
              0, exno(2), exno(1) ];% Such that sigma.n = nomat*sier

    % Dirichlet data
    if order==1
       urr1 = u( [2*no1-1,2*no1,2*no2-1,2*no2], 1 );
       urr2 = u( [2*no1-1,2*no1,2*no2-1,2*no2], 2 );
    elseif order==2
       no3 = bonod(4);
       x3  = nodes(no3,1); y3 = nodes(no3,2);
       urr1 = u( [2*no1-1,2*no1,2*no2-1,2*no2,2*no3-1,2*no3], 1 );
       urr2 = u( [2*no1-1,2*no1,2*no2-1,2*no2,2*no3-1,2*no3], 2 );
    end
    [ Xg, Wg ] = gaussPt1d( Ng );
              
    % Neuman  Data
    fer1 = [f([2*i-1,2*i],1)];
    fer2 = [f([2*i-1,2*i],2)];

    for j=1:Ng
       xg = Xg(j); wg = Wg(j);
       
       % Interpolations
       if order==1
          uer1 = (1-xg)*urr1(1:2) + xg*urr1(3:4); % [ux;uy] on the
          uer2 = (1-xg)*urr2(1:2) + xg*urr2(3:4); % Gauss point
          xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
       elseif order==2
          uer1 = urr1(1:2) + ...
                 xg*(4*urr1(5:6)-3*urr1(1:2)-urr1(3:4)) + ...   % [ux;uy] on the
                 xg^2*(2*urr1(3:4)+2*urr1(1:2)-4*urr1(5:6));    % Gauss point
          uer2 = urr2(1:2) + ...
                 xg*(4*urr2(5:6)-3*urr2(1:2)-urr2(3:4)) + ...  
                 xg^2*(2*urr2(3:4)+2*urr2(1:2)-4*urr2(5:6));
          xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                 xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
       end

        % Test fields
       s1 = [1;0;0]; s2 = [0;1;0]; s3 = [0;0;.5];
       f1 = nomat*s1; f2 = nomat*s2; f3 = nomat*s3;
       v1 = 1/E*[ xgr(1) ; -nu*xgr(2) ];
       v2 = 1/E*[ -nu*xgr(1) ; xgr(2) ];
       v3 = (1+nu)/(2*E)*[ xgr(2) ; xgr(1) ];

       R11 = R11 + len * wg * ( fer1'*v1 - f1'*uer1 ); % increment the integral
       R21 = R21 + len * wg * ( fer1'*v2 - f2'*uer1 );
       R31 = R31 + len * wg * ( fer1'*v3 - f3'*uer1 );
       R12 = R12 + len * wg * ( fer2'*v1 - f1'*uer2 );
       R22 = R22 + len * wg * ( fer2'*v2 - f2'*uer2 );
       R32 = R32 + len * wg * ( fer2'*v3 - f3'*uer2 );
    end
 end

 % Normalize R1
 normR1 = sqrt( 2*(R11^2+R21^2+2*R31^2) - (R11+R21)^2 );
 R1b1 = R11/normR1; R2b1 = R21/normR1; R3b1 = R31/normR1;

 lam11 = (1+R1b1+R2b1)/2; lam21 = -(1-R1b1-R2b1)/2; R1 = [R11,R31;R31,R21];
 [phi1,Lambda1] = eig( [R1b1,R3b1;R3b1,R2b1] );
 dir11 = [ sqrt( abs(Lambda1(1,1)) ) ; sqrt( abs(Lambda1(2,2)) ) ];
 dir21 = [ sign(Lambda1(1,1)) * sqrt( abs(Lambda1(1,1)) ) ;...
           sign(Lambda1(2,2)) * sqrt( abs(Lambda1(2,2)) ) ];
 dir31 = [ sign(Lambda1(1,1)) * sqrt( abs(Lambda1(1,1)) ) ;...
           sqrt( abs(Lambda1(2,2)) ) ];
 dir41 = [ sqrt( abs(Lambda1(1,1)) ) ;...
           sign(Lambda1(2,2)) * sqrt( abs(Lambda1(2,2)) ) ];
          
 dir11 = phi1*dir11; dir11 = dir11/norm(dir11);
 dir21 = phi1*dir21; dir21 = dir21/norm(dir21); % Normal candidates (1)
 dir31 = phi1*dir31; dir31 = dir31/norm(dir31);
 dir41 = phi1*dir41; dir41 = dir41/norm(dir41);

 % Normalize R2
 normR2 = sqrt( 2*(R12^2+R22^2+2*R32^2) - (R12+R22)^2 );
 R1b2 = R12/normR2; R2b2 = R22/normR2; R3b2 = R32/normR2;

 R2 = [R12,R32;R32,R22];
 [phi2,Lambda2] = eig( [R1b2,R3b2;R3b2,R2b2] );
 dir12 = [ sqrt( abs(Lambda2(1,1)) ) ; sqrt( abs(Lambda2(2,2)) ) ];
 dir22 = [ sign(Lambda2(1,1)) * sqrt( abs(Lambda2(1,1)) ) ;...
          sign(Lambda2(2,2)) * sqrt( abs(Lambda2(2,2)) ) ];
 dir12 = phi2*dir12; dir12 = dir12/norm(dir12);
 dir22 = phi2*dir22; dir22 = dir22/norm(dir22); % Normal candidates (2)
 dir32 = [ sign(Lambda2(1,1)) * sqrt( abs(Lambda2(1,1)) ) ;...
           sqrt( abs(Lambda2(2,2)) ) ];
 dir42 = [ sqrt( abs(Lambda2(1,1)) ) ;...
           sign(Lambda2(2,2)) * sqrt( abs(Lambda2(2,2)) ) ];

 % Find the real normal (closest candidates)
 dist11 = norm(dir11-dir12); dist12 = norm(dir11-dir22);
 dist21 = norm(dir21-dir12); dist22 = norm(dir21-dir22);
 dist11m = norm(dir11+dir12); dist12m = norm(dir11+dir22);
 dist21m = norm(dir21+dir12); dist22m = norm(dir21+dir22);
 mindist = min([dist11,dist12,dist21,dist22,dist11m,dist12m,dist21m,dist22m]);

 if dist11 == mindist
    normal = .5*(dir11+dir12);
 elseif dist22 == mindist
    normal = .5*(dir21+dir22);
 elseif dist12 == mindist
    normal = .5*(dir11+dir22);
 elseif dist21 == mindist
    normal = .5*(dir21+dir12);
 elseif dist11m == mindist
    normal = .5*(dir11-dir12);
 elseif dist22m == mindist
    normal = .5*(dir21-dir22);
 elseif dist12m == mindist
    normal = .5*(dir11-dir22);
 elseif dist21m == mindist
    normal = .5*(dir21-dir12);
 end

 % Base matrix
 Q = [ normal(2), normal(1) ; - normal(1), normal(2) ];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Second Task : find the constant

 if order==1
    Ng = 2;
 elseif order==2
    Ng  = 2; 
 end

 norep = Q'*nodes'; K = min(norep(2,:));

 % Norm of [[ut]]
 normT  = sqrt( abs( 2*(R11^2+R21^2+2*R31^2) - 2*(R11+R21)^2 ) );
 normT2 = sqrt( abs( 2*(R12^2+R22^2+2*R32^2) - 2*(R12+R22)^2 ) );

 Rt = 0; Rt2 = 0;
 for i=1:nboun
    bonod = boundary(i,:); exno = extnorm(i,:)';
 
    no1 = bonod(2); no2 = bonod(3);
    x1 = nodes(no1,1); y1 = nodes(no1,2);
    x2 = nodes(no2,1); y2 = nodes(no2,2);
    len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
    
 %   % Compute sigma.n
    nomat = [ exno(1), 0, exno(2) ; ...
              0, exno(2), exno(1) ];% Such that sigma.n = nomat*sier

    % Dirichlet data
    if order==1
       urr1 = u( [2*no1-1,2*no1,2*no2-1,2*no2], 1 );
       urr2 = u( [2*no1-1,2*no1,2*no2-1,2*no2], 2 );
    elseif order==2
       no3 = bonod(4);
       x3  = nodes(no3,1); y3 = nodes(no3,2);
       urr1 = u( [2*no1-1,2*no1,2*no2-1,2*no2,2*no3-1,2*no3], 1 );
       urr2 = u( [2*no1-1,2*no1,2*no2-1,2*no2,2*no3-1,2*no3], 2 );
    end
    [ Xg, Wg ] = gaussPt1d( Ng );
              
    % Neuman  Data
    fer1 = [f([2*i-1,2*i],1)];
    fer2 = [f([2*i-1,2*i],2)];

    for j=1:Ng
       xg = Xg(j); wg = Wg(j);
       
       % Interpolations
       if order==1
          uer1 = (1-xg)*urr1(1:2) + xg*urr1(3:4); % [ux;uy] on the
          uer2 = (1-xg)*urr2(1:2) + xg*urr2(3:4); % Gauss point
          xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
       elseif order==2
          uer1 = urr1(1:2) + ...
                 xg*(4*urr1(5:6)-3*urr1(1:2)-urr1(3:4)) + ...   % [ux;uy] on the
                 xg^2*(2*urr1(3:4)+2*urr1(1:2)-4*urr1(5:6));    % Gauss point
          uer2 = urr2(1:2) + ...
                 xg*(4*urr2(5:6)-3*urr2(1:2)-urr2(3:4)) + ...  
                 xg^2*(2*urr2(3:4)+2*urr2(1:2)-4*urr2(5:6));
          xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                 xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
       end

       % Test fields
       ixigrec = Q'*[xgr(1);xgr(2)]; X = ixigrec(1); Y = ixigrec(2)-K;
       
       sloc = [-X,Y;Y,0];
       st = Q*sloc*Q';
       ft = st*exno;
      
       vloc = [ -X^2/(2*E) + (2+nu)*Y^2/(2*E) ; nu*X*Y/E ];
       vt = Q*vloc;

       Rt  = Rt + len * wg * ( fer1'*vt - ft'*uer1 ); % increment the integral
       Rt2 = Rt2 + len * wg * ( fer2'*vt - ft'*uer2 );
    end
 end

 Cte1 = min( Rt/normT, -Rt/normT ) - K;       % Select the negative one
 Cte2 = min( Rt2/normT2, -Rt2/normT2 ) - K;  %  /!\ The sign depends on the test case

 if cst == 1
    Cte = Cte1;
 elseif cst == 2
    Cte = Cte2;
 else
   error('No of the test-case for the constant computation is not conform');
 end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Last task : identify the displacement gap

 if order==1
    Ng = 2;
 elseif order==2
    Ng  = 2; 
 end

 % Compute L  by the size of the rectangle box
 tangent = [normal(2) ; -normal(1)]; b = max(nodes(:,1)) - min(nodes(:,1));
 a = tangent(2)/tangent(1)*b; L = sqrt(a^2+b^2);
% xmin = min(nodes(:,1)); xmax = max(nodes(:,1));
% ymin = min(nodes(:,2)); ymax = max(nodes(:,2));
% % 4 intersections of the line with the boundaries
% t1 = (xmin-x1)/(x2-x1); t2 = (xmax-x1)/(x2-x1);
% t3 = (ymin-y1)/(y2-y1); t4 = (ymax-y1)/(y2-y1);
% xy11 = [xmin;y1+(y2-y1)*t1];
% xy22 = [xmax;y1+(y2-y1)*t2];
% xy33 = [x1+(x2-x1)*t3;ymin];
% xy44 = [x1+(x2-x1)*t4;ymax];
% xy1234 = [xy11,xy22,xy33,xy44];
% % limit to those that are inside the square
% elim1 = find(xy1234(1,:)>1); elim2 = find(xy1234(2,:)>1);
% elim3 = find(xy1234(1,:)<0); elim4 = find(xy1234(2,:)<0);
% total = setdiff( [1,2,3,4] , union(union(elim1,elim2),union(elim3,elim4)) );
% xyt = xy1234(:,total); % Normally, this one should be of size 2
% xy1 = xyt(:,1); xy2 = xyt(:,2); L = norm(xy1-xy2);
 offset = Q(1,2)/Q(1,1)*Cte;

 % Transform into spacial space
 x1 = offset; x2 = offset+L; xstep = L/100; x = x1:xstep:x2;

 if fourier == 0
    %% Und los geht's !
    % First, build the polynomial test functions.
    coef = zeros(ordp+2, ordp+1);
    coef( 1:2, 1 ) = L*[-(1-nu^2)/(nu*E) ; 0]; % L* because of the homotecy
   
    for k = 1:ordp
       Rhsco = zeros(k+1,1);
       Lhsco = zeros(k+1);
      
       Axx = zeros(k-1,2);
       Axy = zeros(k-1,2);
       Ayy = zeros(k-1,2);
       Bxx = zeros(k-1,2);
       Bxy = zeros(k-1,2);
       Byy = zeros(k-1,2);
      
       azero = -(1-nu^2)/(nu*E*(k+1));
      
       for i=0:floor( (k-1)/2 )
          Lhsco(2*i+1,2*i+1) = (k+1-2*i)*(k-2*i)*E/(1-nu^2);  % a_i
          Lhsco(2*i+1,2*i+2) = (2*i+1)*(k-2*i) * ( nu*E/(1-nu^2) + E/(2*(1+nu)) ); % b_i
          Lhsco(2*i+1,2*i+3) = (2*i+1)*(2*i+2)*E/(2*(1+nu));  % a_{i+1}
       end
      
       for i=1:floor( k/2 )
          Lhsco(2*i,2*i)   = (k-2*i+2)*(k-2*i+1)*E/(2*(1+nu)) ; % b_{i-1}
          Lhsco(2*i,2*i+1) = 2*i*(k+1-2*i)*( E/(2*(1+nu)) + nu*E/(1-nu^2) ) ; % a_i
          Lhsco(2*i,2*i+2) = 2*i*(2*i+1)*E/(1-nu^2) ; %b_i
       end
   
       C = [ eye(2) , zeros(2,k)]; %impose a0 = a0 and b0 = 0 ...
       Lhsco = [ Lhsco ; C ]; % ... in order to find an unique solution
       Rhsco = L*[ Rhsco ; azero ; 0 ];% Alternatively, you can remove ..
         % the L* on Rhsco and Lhs (it's because of the variable change )
      
       Lhsco(size(Lhsco,1)-2,:) = [];  % 'cause square matrix is life
       Rhsco(size(Rhsco,1)-2,:) = [];
      
       coef( 1:k+2, k+1 ) = Lhsco\Rhsco;
    end

    % Place zeros in coef
    coefa = coef;
    coefb = coef;
    for i=1:size(coefa,1)
       if mod(i,2) == 0
          coefa(i,:) = 0;
       end
       if mod(i,2) == 1
          coefb(i,:) = 0;
       end
    end
   
    %% Compute the RG
    Rhs1  = zeros(ordp+1,1); Rhs2  = zeros(ordp+1,1);
    vpa  = zeros(2*nnodes, 1);

    for k=0:ordp
       Rhs1(k+1) = 0; Rhs2(k+1) = 0;
       for i=1:nboun % boundary1 and boundary 2 are supposed to be the same
          bonod = boundary(i,:); exno = extnorm(i,:)';
    
          no1 = bonod(2); no2 = bonod(3);
          x1 = nodes(no1,1); y1 = nodes(no1,2);
          x2 = nodes(no2,1); y2 = nodes(no2,2);
          len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
    
       %   % Compute sigma.n
          nomat = [ exno(1), 0, exno(2) ; ...
                    0, exno(2), exno(1) ];% Such that sigma.n = nomat*sier
   
          % Dirichlet data
          if order==1
             urr1 = u( [2*no1-1,2*no1,2*no2-1,2*no2], 1 );
             urr2 = u( [2*no1-1,2*no1,2*no2-1,2*no2], 2 );
          elseif order==2
             no3 = bonod(4);
             x3  = nodes(no3,1); y3 = nodes(no3,2);
             urr1 = u( [2*no1-1,2*no1,2*no2-1,2*no2,2*no3-1,2*no3], 1 );
             urr2 = u( [2*no1-1,2*no1,2*no2-1,2*no2,2*no3-1,2*no3], 2 );
          end
          [ Xg, Wg ] = gaussPt1d( Ng );
              
          % Neuman  Data
          fer1 = [f([2*i-1,2*i],1)];
          fer2 = [f([2*i-1,2*i],2)];
   
          for j=1:Ng
             xg = Xg(j); wg = Wg(j);
          
             % Interpolations
             if order==1
                uer1 = (1-xg)*urr1(1:2) + xg*urr1(3:4); % [ux;uy] on the
                uer2 = (1-xg)*urr2(1:2) + xg*urr2(3:4); % Gauss point
                xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
             elseif order==2
                uer1 = urr1(1:2) + ...
                       xg*(4*urr1(5:6)-3*urr1(1:2)-urr1(3:4)) + ...   % [ux;uy] on the
                       xg^2*(2*urr1(3:4)+2*urr1(1:2)-4*urr1(5:6));    % Gauss point
                uer2 = urr2(1:2) + ...
                       xg*(4*urr2(5:6)-3*urr2(1:2)-urr2(3:4)) + ...  
                       xg^2*(2*urr2(3:4)+2*urr2(1:2)-4*urr2(5:6));
                xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                       xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
             end
            
             ixigrec = Q'*[xgr(1);xgr(2)];
             X = (ixigrec(1)-offset)/L-.5; Y = (ixigrec(2)+Cte)/L;
            
             % Build X^k+1-j*Y^j
             GROX = zeros(ordp+2,1);
             for j = 0:k+1
                GROX(j+1) = X^(k+1-j)*Y^j;
             end
            
             vloc = [ coefa(:,k+1)'*GROX ; coefb(:,k+1)'*GROX ];
             vpa = Q*vloc;
          
             sloc11 = 0; sloc12 = 0; sloc22 = 0;
             for j=0:floor(k/2)
                sloc11 = sloc11 + ( E/(1-nu^2)*(k+1-2*j)*coefa(2*j+1,k+1) + ...
                         nu*E/(1-nu^2)*(2*j+1)*coefb(2*j+2,k+1) )*X^(k-2*j)*Y^(2*j);
                sloc22 = sloc22 + ( nu*E/(1-nu^2)*(k+1-2*j)*coefa(2*j+1,k+1) + ...
                         E/(1-nu^2)*(2*j+1)*coefb(2*j+2,k+1) )*X^(k-2*j)*Y^(2*j);
             end
             for j=1:floor((k+1)/2)
                sloc12 = sloc12 + ( E/(2*(1+nu))*(2*j)*coefa(2*j+1,k+1)* ...
                                    X^(k+1-2*j)*Y^(2*j-1) );
             end
             for j=0:floor((k-1)/2)
                sloc12 = sloc12 + ( E/(2*(1+nu))*(k-2*j)*coefb(2*j+2,k+1)* ...
                                    X^(k-1-2*j)*Y^(2*j+1) );
             end
   
             sloc = 1/L*[sloc11,sloc12;sloc12,sloc22];
             spa = Q*sloc*Q';
             fpa = spa*exno;
   
             Rhs1(k+1) = Rhs1(k+1) + len * wg * ( fer1'*vpa - fpa'*uer1 );
             Rhs2(k+1) = Rhs2(k+1) + len * wg * ( fer2'*vpa - fpa'*uer2 );
          end
       end
      
    end

    L1 = -.5; L2 = .5;
    for i=0:ordp
       for j=0:ordp
          ord = i+j+1;
          Lhs(i+1,j+1) = L*(L2^ord - L1^ord)/ord;
       end
    end
   
    % Invert the operator
    McCoef = Lhs\[Rhs1,Rhs2];

    nbase = size(McCoef,1);
   
    ug = zeros(size(x,2),2);
    for i=1:nbase
       ug(:,1) = ug(:,1) + McCoef(i,1)*( (x'-offset)./L-.5).^(i-1);
       ug(:,2) = ug(:,2) + McCoef(i,2)*( (x'-offset)./L-.5).^(i-1);
    end

 else % Use Fourier
    nbase = ordp;

    Rp     = zeros(nbase+1,1);
    Rm     = zeros(nbase+1,1);
    Rpe    = zeros(nbase+1,1);
    Rme    = zeros(nbase+1,1);
   
    lambda = zeros(nbase+1,1);
    fourn  = zeros(nbase+1,1);
    fournm = zeros(nbase+1,1);
    akan   = zeros(nbase,1);
    bkan   = zeros(nbase,1);
   
    for kp=2:nbase+1
       k = kp-1;
       vp = zeros(2*nnodes,1);
       vm = zeros(2*nnodes,1);
       lambda(kp) = 2*k*pi/L;
       for sx = [1,-1]  % sx=-1 is not really used, but Debug stuff
          lambdae = sx*lambda(kp);
   
          if sx == 1 Rp(kp) = 0; else Rpm = 0; end
          for i=1:nboun % boundary1 and boundary 2 are supposed to be the same
             bonod = boundary(i,:); exno = extnorm(i,:)';
    
             no1 = bonod(2); no2 = bonod(3);
             x1 = nodes(no1,1); y1 = nodes(no1,2);
             x2 = nodes(no2,1); y2 = nodes(no2,2);
             len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
    
          %   % Compute sigma.n
             nomat = [ exno(1), 0, exno(2) ; ...
                       0, exno(2), exno(1) ];% Such that sigma.n = nomat*sier
   
             % Dirichlet data
             if order==1
                urr1 = u( [2*no1-1,2*no1,2*no2-1,2*no2], 1 );
                urr2 = u( [2*no1-1,2*no1,2*no2-1,2*no2], 2 );
             elseif order==2
                no3 = bonod(4);
                x3  = nodes(no3,1); y3 = nodes(no3,2);
                urr1 = u( [2*no1-1,2*no1,2*no2-1,2*no2,2*no3-1,2*no3], 1 );
                urr2 = u( [2*no1-1,2*no1,2*no2-1,2*no2,2*no3-1,2*no3], 2 );
             end
             [ Xg, Wg ] = gaussPt1d( Ng );
              
             % Neuman  Data
             fer1 = [f([2*i-1,2*i],1)];
             fer2 = [f([2*i-1,2*i],2)];
   
             for j=1:Ng
                xg = Xg(j); wg = Wg(j);
          
                % Interpolations
                if order==1
                   uer1 = (1-xg)*urr1(1:2) + xg*urr1(3:4); % [ux;uy] on the
                   uer2 = (1-xg)*urr2(1:2) + xg*urr2(3:4); % Gauss point
                   xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
                elseif order==2
                   uer1 = urr1(1:2) + ...
                          xg*(4*urr1(5:6)-3*urr1(1:2)-urr1(3:4)) + ...   % [ux;uy] on the
                          xg^2*(2*urr1(3:4)+2*urr1(1:2)-4*urr1(5:6));    % Gauss point
                   uer2 = urr2(1:2) + ...
                          xg*(4*urr2(5:6)-3*urr2(1:2)-urr2(3:4)) + ...  
                          xg^2*(2*urr2(3:4)+2*urr2(1:2)-4*urr2(5:6));
                   xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                          xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
                end
         
                % Test fields
                ixigrec = Q'*[xgr(1);xgr(2)]; X = ixigrec(1); Y = ixigrec(2)+Cte;
                
                sloc11 = E/(1+nu)*lambdae^2*exp(-I*lambdae*X)* ...
                                ( exp(lambda(kp)*Y)+exp(-lambda(kp)*Y) );
                sloc12 = E/(1+nu)*I*lambdae^2*exp(-I*lambdae*X)* ...
                                   ( exp(lambda(kp)*Y)-exp(-lambda(kp)*Y) );

                sloc = [-sloc11,-sloc12;-sloc12,sloc11];
                sp = Q*sloc*Q';
                fp = sp*exno;
               
                v1 = -I*lambda(kp)*exp(-I*lambdae*X)* ...
                                ( exp(lambda(kp)*Y)+exp(-lambda(kp)*Y) );
                v2 = lambda(kp)*exp(-I*lambdae*X)* ...
                                ( exp(lambda(kp)*Y)-exp(-lambda(kp)*Y) );
                vloc = [ v1 ; v2 ];
                vp = Q*vloc;

                if sx == 1
                   Rp(kp) = Rp(kp) + len * wg * ( fer1'*vp - fp'*uer1 );
                   fourn(kp) = - (1+nu)/(2*E*L*lambda(kp)^2)*Rp(kp);
                   akan(k) = 2*real(fourn(kp));
                   bkan(k) = 2*imag(fourn(kp));
                else
                   Rpm = Rpm + len * wg * ( fer1'*vp - fp'*uer1 );
                   fournm(kp) = -(1+nu)/(2*E*L*lambda(kp)^2)*Rpm;
                end  
             end
          end
               
       end
    end
   
    % The constant term
    fourn(1) = -(R11+R21)/L;
   
    % Reconstruct the normal displacement
    ug = fourn(1) + sum ( [0*x ;...  % Hack in case there is 1 element only (to let the sum work)
            akan.*cos(lambda(2:end)*x) + bkan.*sin(lambda(2:end)*x) ] );
  end

end
