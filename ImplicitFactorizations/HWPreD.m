function z = HWPreD(K,M,n,gamma,alpha,Dvort,Dn,r)
%block diagonal preconitioner
z=zeros(3*n,1);
z(1:n) = (M - gamma*Dvort*K)\r(1:n);

z(n+1:2*n) = - (M\r(n+1:2*n))/(alpha*gamma);


z(2*n+1:3*n) = r(2*n+1:3*n);