
% 2D Quadrilateral Finite Element Methods Using Discontinuous Pressures
clc
clear

nu=1;
f=@(x,y) sin(x)+sin(y);

x_max=1;
x_min=0;
y_max=1;
y_min=0;
x_n=3;
y_n=2;
node_number=(x_n+1)*(y_n+1);
element_number=x_n*y_n;
ubasis_function_number=x_n*(y_n+1)+y_n*(x_n+1)+node_number*2;
pbasis_function_number=element_number;
x=x_min:(x_max-x_min)/x_n:x_max;
y=y_min:(y_max-y_min)/y_n:y_max;
[X,Y]=meshgrid(x,y);
X=reshape(X',[],1);
Y=reshape(Y',[],1);
node_coordinate=[X,Y];
[m,n]=size(X);
for i=1:m
    if(node_coordinate(i,2)==0||node_coordinate(i,1)==0||node_coordinate(i,1)==1)
        node_coordinate(i,3)=1;  % 1 = boundary 0= interior
         node_coordinate(i,4)=0;  % Dirichlet boundary
    elseif(node_coordinate(i,2)==1)
         node_coordinate(i,3)=1;  % 1 = boundary 0= interior
         node_coordinate(i,4)=1;  % Dirichlet boundary
    else
        node_coordinate(i,3)=0;  % 1 = boundary 0= interior
         node_coordinate(i,4)=0;  % Dirichlet boundary
    end
end

element_coordinate11=(1:x_n)';
edge_coordinate11=(1:x_n)';
edge_coordinate1=[];
element_coordinate1=[];
for i=1:y_n
    element_coordinate1=[element_coordinate1;element_coordinate11];
    element_coordinate11=element_coordinate11+x_n+1;
    edge_coordinate1=[edge_coordinate1;edge_coordinate11];
    edge_coordinate11=edge_coordinate11+2*x_n+1;
end
element_coordinate=[element_coordinate1,element_coordinate1+1,element_coordinate1+x_n+2,element_coordinate1+x_n+1,edge_coordinate1,edge_coordinate1+x_n+1,edge_coordinate1+2*x_n+1,edge_coordinate1+x_n];



% phi1=@(x,y)[1-x-y+x.*y-3*y.*(1-x).*(1-y);0];
% phi2=@(x,y)[x-x.*y-3*x.*y.*(1-y);0];
% phi3=@(x,y)[x.*y-3*x.*y.*(1-y);0];
% phi4=@(x,y)[y-x.*y-3*y.*(1-x).*(1-y);0];
% phi5=@(x,y)[-6*y.*(1-x).*(1-y);0];
% phi6=@(x,y)[6*x.*y.*(1-y);0];
% 
% phi7=@(x,y)[0;1-x-y+x.*y-3*x.*(1-x).*(1-y)];
% phi8=@(x,y)[0;x-x.*y-3*x.*(1-x).*(1-y)];
% phi9=@(x,y)[0;x.*y-3*x.*y.*(1-x)];
% phi10=@(x,y)[0;y-x.*y-3*x.*y*(1-x)];
% phi11=@(x,y)[0;-6*x.*(1-x).*(1-y)];
% phi12=@(x,y)[0;6*x.*y.*(1-x)];
syms x y
coefmat=coef();
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

% local matrix A 
% for i=1:12
%     for j=1:i
%         A11=int(matlabFunction(phix(1,i).*phix(1,j),'vars',[x y]),x,-1,1);
%         A1=int(matlabFunction(A11,'var',[y]),y,-1,1);
%         A22=int(matlabFunction(phix(2,i).*phix(2,j),'vars',[x y]),x,-1,1);
%         A2=int(matlabFunction(A22,'var',[y]),y,-1,1);
%         A33=int(matlabFunction(phiy(1,i).*phiy(1,j),'vars',[x y]),x,-1,1);
%         A3=int(matlabFunction(A33,'var',[y]),y,-1,1);
%         A44=int(matlabFunction(phiy(2,i).*phiy(2,j),'vars',[x y]),x,-1,1);
%         A4=int(matlabFunction(A44,'var',[y]),y,-1,1);
%         A(i,j)=double(A1+A2+A3+A4);
%         A(j,i)=A(i,j);
%     
%     end
% end
load A11reference.mat;
load A12reference.mat;
load A22reference.mat;
load B11reference.mat;
load B12reference.mat;
load B22reference.mat;

A=zeros(ubasis_function_number);

for i=1:element_number 
    local_coordinate1=node_coordinate(element_coordinate(i,1),1:2);
    local_coordinate2=node_coordinate(element_coordinate(i,2),1:2);
    local_coordinate3=node_coordinate(element_coordinate(i,3),1:2);
    local_coordinate4=node_coordinate(element_coordinate(i,4),1:2);
    J=[(local_coordinate2-local_coordinate1).'/2,(local_coordinate3-local_coordinate2).'/2];
    detJ=det(J);
    invJ=inv(J);
    invJT=invJ.';
    C=invJ*invJT;
    
    Alocal=nu*detJ*(C(1,1)*(A11reference+B11reference)+C(1,2)*(A12reference+A12reference.'+B12reference+B12reference.')+C(2,2)*(A22reference+B22reference));

    for j=1:12
        for k=1:12
            if j<=8
                jglobal2=mod(j-1,2)+1;
                jglobal1=element_coordinate(i,floor((j-1)/2)+1);
                jglobal=2*jglobal1+jglobal2-2;
               
            else
                jglobal1=j-8;
                jglobal=2*node_number+element_coordinate(i,jglobal1+4);
                
            end
            if k<=8
                kglobal2=mod(k-1,2)+1;
                kglobal1=element_coordinate(i,floor((k-1)/2)+1);
                kglobal=2*kglobal1+kglobal2-2;
               
            else
                kglobal1=k-8;
                kglobal=2*node_number+element_coordinate(i,kglobal1+4);
              
            end
            A(jglobal,kglobal)=A(jglobal,kglobal)+Alocal(j,k);
        end
    end
end