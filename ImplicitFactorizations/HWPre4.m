function z = HWPre4(P,K,M,n,gamma,alpha,Dvort,Dn,r)
% Preconditioner P4 for the Hasegawa-Wakatani problem

Dw=Dvort;
z1=zeros(3*n,1);
z1(1:n) = r(2*n+1:3*n);
z1(2*n+1:3*n) = r(1:n) + z1(1:n);
z1(n+1:2*n) = r(n+1:2*n) - (gamma*Dw*Dn/alpha)*(K*(M\(K*(M\z1(1:n))))) - (1/(alpha*gamma))*(M + gamma*(alpha*M-Dn*K))*(M\z1(2*n+1:3*n)) ;

z2=zeros(3*n,1);
z2(2*n+1:3*n) = -M\z1(1:n);
temp = z1(n+1:2*n) - ((Dw*(1+alpha*gamma))/alpha)*K*z2(2*n+1:3*n);
z2(n+1:2*n) = -(alpha*gamma)*((((1+alpha*gamma)*M-gamma*Dn*K)/M)*K)\temp;
z2(1:n) = (M\(z1(2*n+1:3*n) - K*z2(n+1:2*n) +gamma*Dw*(K*z2(2*n+1:3*n))))/(alpha*gamma);

z=zeros(3*n,1);
z(n+1:2*n) = z2(n+1:2*n);
z(1:n) = z2(2*n+1:3*n) +M\(K*z(n+1:2*n));
z(2*n+1:3*n) = z2(1:n) + (Dw/alpha)*(M\(K*(M\(K*z(n+1:2*n)))));