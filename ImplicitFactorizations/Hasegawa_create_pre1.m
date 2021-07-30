function P= Hasegawa_create_pre1(K,M,gamma)

n = size(K,2);
alpha=1.0;
%kappa = 0.5;
Dvort = 1.e-4;
Dn = 1.e-4;
B = [-alpha*M + Dn*K, 0*M, alpha*M; -alpha*M, Dvort*K, 0*M; 0*M, 0*M, -0*K ];
C = sparse(3*n,3*n);
C(1:n,1:n) = M;
C(n+1:2*n,n+1:2*n) = M;
C(2*n+1:3*n,n+1:2*n) = -M;
C(2*n+1:3*n,2*n+1:3*n) = K;

A = C-gamma*B;

P = A([n+1:2*n,1:n,2*n+1:3*n],[n+1:3*n, 1:n ]);
P(1:2*n,1:2*n)=speye(2*n);
