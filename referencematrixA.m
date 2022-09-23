function [A11reference,A12reference,A22reference,B11reference,B12reference,B22reference]=referencematrixA(coefmat)
syms x y
phi=[];
phi1var=[1 ;x; y; x.*y ;y.*(1-x).*(1-y); x.*y.*(1-y)];
phi2var=[1 ;x; y; x.*y ;x.*(1-x).*(1-y); x.*y.*(1-x)];
for i=1:12
    phi1=coefmat(i,1:6)*phi1var;
    phi2=coefmat(i,7:12)*phi2var;
    phicol=[phi1;phi2];
    phi=[phi,phicol];
end

phix=diff(phi,x);
phiy=diff(phi,y);

gaussweight=[5/9,8/9,5/9];
gausspoint=[-sqrt(3/5),0,sqrt(3/5)];

% reference matrix A 
for i=1:12
    for j=1:12
        A11=int(matlabFunction(phix(1,i).*phix(1,j),'vars',[x y]),x,-1,1);
        A1=int(matlabFunction(A11,'var',[y]),y,-1,1);
        A11reference(i,j)=double(A1);
%         
      

        A22=int(matlabFunction(phiy(1,i).*phiy(1,j),'vars',[x y]),x,-1,1);
        A2=int(matlabFunction(A22,'var',[y]),y,-1,1);
        A22reference(i,j)=double(A2);
  
          A12=int(matlabFunction(phix(1,i).*phiy(1,j),'vars',[x y]),x,-1,1);
        A3=int(matlabFunction(A12,'var',[y]),y,-1,1);
        A12reference(i,j)=double(A3);


        B11=int(matlabFunction(phix(2,i).*phix(2,j),'vars',[x y]),x,-1,1);
       B1=int(matlabFunction(B11,'var',[y]),y,-1,1);
         B11reference(i,j)=double(B1);
         
        B22=int(matlabFunction(phiy(2,i).*phiy(2,j),'vars',[x y]),x,-1,1);
        B2=int(matlabFunction(B22,'var',[y]),y,-1,1);
        B22reference(i,j)=double(B2);


        B12=int(matlabFunction(phix(2,i).*phiy(2,j),'vars',[x y]),x,-1,1);
        B3=int(matlabFunction(B12,'var',[y]),y,-1,1);
        B12reference(i,j)=double(B3);

    end
end

save A11reference.mat A11reference;
save A12reference.mat A12reference;
save A22reference.mat A22reference;
save B11reference.mat B11reference;
save B12reference.mat B12reference;
save B22reference.mat B22reference;
end