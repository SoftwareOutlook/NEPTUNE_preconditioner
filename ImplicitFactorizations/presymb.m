% Generates non-symmetric implicit-factorization constraint preconditioner
% families
format compact
Ll = [sym('l11'),sym('l12'), sym('l13'),sym('l21'),sym('l22'),...
    sym('l23'),sym('l31'),sym('l32'),sym('l33') ];
Rr = [sym('r11'),sym('r12'), sym('r13'),sym('r21'),sym('r22'),...
    sym('r23'),sym('r31'),sym('r32'),sym('r33') ]
Mm = [sym('m11'),sym('m12'), sym('m13'),sym('m21'),sym('m22'),...
    sym('m23'),sym('m31'),sym('m32'),sym('m33') ]
total = 0;

for i=1:5
    for j=1:3
        for k=1:5
            L = [Ll(1:3);Ll(4:6);Ll(7:9)];
            R = [Rr(1:3);Rr(4:6);Rr(7:9)];
            M = [Mm(1:3);Mm(4:6);Mm(7:9)];
            switch i
                case 1
                    L(1,1:2)=0; L(2,1)=0;
                case 2
                    L(1,1:2)=0; L(3,2)=0;
                case 3
                    L(1,2)=0; L(3,1:2)=0;
                case 4
                    L(2,1)=0; L(3,1:2)=0;    
                case 5
                    L(1,2)=0; L(3,2:3)=0;
            end
            
            switch j
                case 1
                    M(3,2:3)=0; M(2,3)=0;
                case 2
                    M(1,1:2)=0; M(2,1)=0;
                case 3
                    M(1,2)=0; M(2,1)=0; M(2,3)=0; M(3,2)=0;   
            end
            
            switch k
                case 1
                    R(1:2,1)=0; R(1,2)=0;
                case 2
                    R(1:2,1)=0; R(2,3)=0;
                case 3
                    R(2,1)=0; R(2:3,3)=0;
                case 4
                    R(2,1)=0; R(1:2,3)=0;
                case 5
                    R(1,2)=0; R(1:2,3)=0;
            end
            p = 1;
            F = L*M*R;
            if ((F(1,3)==0) | (F(2,3)==0) | (F(3,1)==0) | (F(3,2)==0))
                p=0;
            end
            
            if (p==1)
                total = total+1;
               % [i,j,k]
               factor = total
                struct=[L,M,R]
                F
            end
        end
    end
end
total