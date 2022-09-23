function coefmat=coef()
% syms a1 a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 c4 x1 x2
% 
% x3=1-x1;
% x4=1-x2;
% f1=a1+a2.*x1+a3.*x2+a4.*x1.*x2-c1.*x2.*x3.*x4+c3.*x1.*x2.*x4;
% f2=b1+b2.*x1+b3.*x2+b4.*x1.*x2-c2.*x1.*x3.*x4+c4.*x1.*x2.*x3;
% a11=subs(f1,[x1,x2],[-1,-1])
% a12=subs(f2,[x1,x2],[-1,-1])
% a21=subs(f1,[x1,x2],[1,-1])
% a22=subs(f2,[x1,x2],[1,-1])
% a31=subs(f1,[x1,x2],[1,1])
% a32=subs(f2,[x1,x2],[1,1])
% a41=subs(f1,[x1,x2],[-1,1])
% a42=subs(f2,[x1,x2],[-1,1])
% 
% int1p=subs(f2,x2,-1);
% int11=int(int1p,x1,-1,1)
% int2p=subs(f1,x1,1);
% int21=int(int2p,x2,-1,1)
% int3p=subs(f2,x2,1);
% int31=int(int3p,x1,-1,1)
% int4p=subs(f1,x1,-1);
% int41=int(int4p,x2,-1,1)
A=[1 ,-1, -1 ,1 ,4 ,2 ,0, 0, 0, 0, 0 ,0;
    0 ,0, 0, 0 ,0 ,0, 1 ,-1 ,-1, 1 ,4 ,2;
    1 ,1 ,-1, -1, 0, -2, 0, 0, 0, 0, 0, 0 ;
     0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,-1 ,-1 ,0 ,0 ;
     1 ,1 ,1 ,1 ,0, 0 ,0, 0, 0, 0, 0, 0 ;
     0 ,0 ,0 ,0, 0 ,0 ,1 ,1 ,1 ,1 ,0 ,0;
    1 ,-1 ,1 -1 0 ,0 ,0 ,0 ,0, 0 ,0 ,0;
    0 ,0 ,0 ,0, 0,  0, 1, -1 ,1 ,-1, 0, -2;
    0 ,0 ,0 ,0 ,0, 0, 2, 0, -2,  0, 4/3, 2/3;
    2 ,2 0 ,0 ,0 ,-2/3 ,0 ,0 ,0 ,0 ,0 ,0 ;
    0 ,0, 0, 0, 0 ,0 ,2 ,0, 2 ,0 ,0,-2/3;
    2 ,-2 ,0 ,0, 4/3 ,2/3, 0 ,0 ,0 ,0 ,0 ,0];
coefmat=[];

for i=1:12
    b=zeros(12,1);
    b(i)=1;
    x=A\b;
    coefmat=[coefmat;x'];
end
