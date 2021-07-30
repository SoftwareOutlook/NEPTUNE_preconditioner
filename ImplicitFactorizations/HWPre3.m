function z = HWPre3(AP,K,M,n,gamma,alpha,Dvort,Dn,r)
% Preconditioner P3 for the Hasegawa-Wakatani problem
 z=zeros(3*n,1);
 z(2*n+1:3*n) = (M\r(1:n))/(alpha*gamma);

 rtt = r(n+1:2*n) -(M-gamma*(Dn*K - alpha*M))*z(2*n+1:3*n);
 z(n+1:2*n) = -(M\rtt)/(alpha*gamma); %G22=-gamma*alpha*M
% 
 rt = r(2*n+1:3*n) - K*z(n+1:2*n);
 z(1:n) = -M\rt;

