function [B11reference,B12reference,B21reference,B22reference]=referencematrixB(phi,phix,phiy)

syms x y;

for i=1:12
        B11=int(matlabFunction(phix(1,i),'vars',[x y]),x,-1,1);
        B1=int(matlabFunction(B11,'var',[y]),y,-1,1);
        B11reference(i,1)=-double(B1);

         B12=int(matlabFunction(phiy(1,i),'vars',[x y]),x,-1,1);
        B2=int(matlabFunction(B12,'var',[y]),y,-1,1);
        B12reference(i,1)=-double(B2);

         B21=int(matlabFunction(phix(2,i),'vars',[x y]),x,-1,1);
        B3=int(matlabFunction(B21,'var',[y]),y,-1,1);
        B21reference(i,1)=-double(B3);

         B22=int(matlabFunction(phiy(2,i),'vars',[x y]),x,-1,1);
        B4=int(matlabFunction(B22,'var',[y]),y,-1,1);
        B22reference(i,1)=-double(B4);
end


save B11reference.mat B11reference;
save B12reference.mat B12reference;
save B21reference.mat B21reference;
save B22reference.mat B22reference;