function z = HWPre2(P,K,M,n,gamma,alpha,Dvort,Dn,r)
% Preconditioner P2 for the Hasegawa-Wakatani problem

 z=zeros(3*n,1);
 z(2*n+1:3*n) = (M\r(1:n))/(alpha*gamma);
  rtt = r(n+1:2*n) -(M-gamma*(Dn*K - alpha*M))*z(2*n+1:3*n);
  z(n+1:2*n) = rtt;%G22=I
  rt = r(2*n+1:3*n) - K*z(n+1:2*n);
  z(1:n) = -M\rt;
