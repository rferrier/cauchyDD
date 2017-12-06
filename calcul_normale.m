clear all

eps=10^(-4); % perturbation

errC = 0;
errP = 0;

for i=1:1
   % fabrication des données
   U1 = rand(3,1);
   U2 = rand(3,1);
   N = rand(3,1); N = N / (norm(N)*sign(N(1)));

   rand1=rand(3,3);
   rand2=rand(3,3);
   
   R1 = U1*N' + N*U1' + (rand1)*eps; 
   R2 = U2*N' + N*U2' + (rand2)*eps;
   
   %%%%%%%%%%%%%%%%%%%%%%
   % Methode classique
   %%%%%%%%%%%%%%%%%%%%%%
   
   [A1,S1,B1] = svd(R1) ; 
   [A2,S2,B2] = svd(R2) ; 
   Nclass = cross(A1(:,3),A2(:,3));
   Nclass = Nclass * sign(Nclass(1))/norm(Nclass);
   
   %%%%%%%%%%%%%%%%%%%%%%%%
   % Peter
   %%%%%%%%%%%%%%%%%%%%%%%%
   alpha1 = sqrt((norm(R1,'fro')^2 - trace(R1)^2/2)/2);
   alpha2 = sqrt((norm(R2,'fro')^2 - trace(R2)^2/2)/2);
   
   %Ut1=U1/alpha1;
   %Ut2=U2/alpha2; % 
   Rt1 = R1/alpha1; % Rtilde
   Rt2 = R2/alpha2;
   
   beta1 = trace(Rt1)/2;
   beta2 = trace(Rt2)/2;
   % beta = N'*Ut;
   
   Rc1 = Rt1^2-trace(Rt1)*Rt1/2; % Rchapeau
   Rc2 = Rt2^2-trace(Rt2)*Rt2/2;
   % Rc1 = Ut1*Ut1' + N*N' 
   % Rc2 = Ut2*Ut2' + N*N' 
   % trace(Rc)=2 ;
   
   Q = Rc1 * Rc2 + Rc2*Rc1 - (Rt1*trace(Rt1)/2) - (Rt2*trace(Rt2)/2) ;
   %Q = (Ut1*Ut2'+Ut2*Ut1')*(Ut1'*Ut2)+2*N*N'
   
   gamma = sqrt((trace(Q)-2)/2); % gamma = Ut1'*Ut2;
   
   A = Q - Rc1 * (beta2*gamma/beta1) - Rc2 * (beta1*gamma/beta2) ;
   
   Delta = A - (2-((beta1/beta2)+(beta2/beta1))*gamma)*eye(3,3);
   
   [U,S,V] = svd(Delta) ; % N
   NN = U(:,3);
   NN = NN*sign(NN(1));
   
   errC = errC + norm(Nclass-N);
   errP = errP + norm(NN-N);
end


disp(['erreur Peter ',num2str(errP)]) % quasi 0
disp(['erreur classique ',num2str(errC)]) 
