
% 2D Quadrilateral Finite Element Methods Using Discontinuous Pressures
clc
clear

nu=1;
syms xi eta 
%f1=1;
%f2=1;
% g1= [sin(pi*xi),0];
% % g2=[1,0];
% g2=[sin(pi*xi),-pi*cos(pi*xi)];
% g3=[sin(pi*eta),pi*eta];
% g4=[sin(pi*eta),-pi*eta];

x_max=1;
x_min=0;
y_max=1;
y_min=0;
x_n=3;
y_n=4;
localbasisfunctionnumber=12;
node_number=(x_n+1)*(y_n+1);
element_number=x_n*y_n;
edge_number=x_n*(y_n+1)+y_n*(x_n+1);
ubasis_function_number=edge_number+node_number*2;
pbasis_function_number=element_number;
deltax=(x_max-x_min)/x_n;
deltay=(y_max-y_min)/y_n;
x=x_min:deltax:x_max;
y=y_min:deltay:y_max;
[X,Y]=meshgrid(x,y);
X=reshape(X',[],1);
Y=reshape(Y',[],1);
node_coordinate=[X,Y];
[m,n]=size(X);

for i=1:m
    if(node_coordinate(i,2)==0||node_coordinate(i,1)==0||node_coordinate(i,2)==1||node_coordinate(i,1)==1)
        node_coordinate(i,3)=1;  % 1 = boundary 0= interior
        node_coordinate(i,4)=function_u(node_coordinate(i,1),node_coordinate(i,2));
        node_coordinate(i,5)=function_v(node_coordinate(i,1),node_coordinate(i,2));
    else
        node_coordinate(i,3)=0;  % 1 = boundary 0= interior
        node_coordinate(i,4)=0;
        node_coordinate(i,5)=0;
    end
end
% for i=1:m
%     if(node_coordinate(i,2)==0)
%         node_coordinate(i,3)=1;  % 1 = boundary 0= interior
%         nodevalue=double(subs(g1,xi,node_coordinate(i,1)));
%          node_coordinate(i,4)=nodevalue(1,1);% Dirichlet boundary
%          node_coordinate(i,5)=nodevalue(1,2);
%     elseif(node_coordinate(i,1)==0)
%         node_coordinate(i,3)=1; 
%         nodevalue=double(subs(g4,eta,node_coordinate(i,2)));
%          node_coordinate(i,4)=nodevalue(1,1);% Dirichlet boundary
%          node_coordinate(i,5)=nodevalue(1,2);
%     elseif(node_coordinate(i,2)==1)
%          node_coordinate(i,3)=1;  % 1 = boundary 0= interior
%          nodevalue=double(subs(g2,xi,node_coordinate(i,1)));
%          node_coordinate(i,4)=nodevalue(1,1);% Dirichlet boundary
%          node_coordinate(i,5)=nodevalue(1,2);
%     elseif(node_coordinate(i,1)==1)
%         node_coordinate(i,3)=1;  % 1 = boundary 0= interior
%          nodevalue=double(subs(g3,eta,node_coordinate(i,2)));
%          node_coordinate(i,4)=nodevalue(1,1);% Dirichlet boundary
%          node_coordinate(i,5)=nodevalue(1,2);
%     else
%         node_coordinate(i,3)=0;  % 1 = boundary 0= interior
%          node_coordinate(i,4)=0;  % Dirichlet boundary
%          node_coordinate(i,5)=0; 
%     end
% end
% % 
% for i=1:m
%     if(node_coordinate(i,2)==0||node_coordinate(i,1)==0||node_coordinate(i,1)==1)
%         node_coordinate(i,3)=1;  % 1 = boundary 0= interior
%          node_coordinate(i,4)=g1(1,1);  % Dirichlet boundary
%          node_coordinate(i,5)=g1(1,2); 
%     elseif(node_coordinate(i,2)==1)
%          node_coordinate(i,3)=1;  % 1 = boundary 0= interior
%          node_coordinate(i,4)=g2(1,1);% Dirichlet boundary
%          node_coordinate(i,5)=g2(1,2);
%     else
%         node_coordinate(i,3)=0;  % 1 = boundary 0= interior
%          node_coordinate(i,4)=0;  % Dirichlet boundary
%          node_coordinate(i,5)=0; 
%     end
% end

% for i=1:m
%     if (node_coordinate(i,2)==1)
%          node_coordinate(i,3)=1;  % 1 = boundary 0= interior
%          node_coordinate(i,4)=g2(1,1);% Dirichlet boundary
%          node_coordinate(i,5)=g2(1,2);
%        
%     elseif (node_coordinate(i,2)==0||node_coordinate(i,1)==0||node_coordinate(i,1)==1)
%          node_coordinate(i,3)=1;  % 1 = boundary 0= interior
%          node_coordinate(i,4)=g1(1,1);  % Dirichlet boundary
%          node_coordinate(i,5)=g1(1,2); 
%     else
%         node_coordinate(i,3)=0;  % 1 = boundary 0= interior
%          node_coordinate(i,4)=0;  % Dirichlet boundary
%          node_coordinate(i,5)=0; 
%     end
% end




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



edge_dirichlet=zeros(edge_number,2);

g1=[function_u(xi,0),function_v(xi,0)];
g2=[function_u(xi,1),function_v(xi,1)];
g3=[function_u(1,eta),function_v(1,eta)];
g4=[function_u(0,eta),function_v(0,eta)];
for i=1:x_n
    intg1=int(g1*[0;-1],xi,(i-1)/x_n,i/x_n);
    intg2=int(g2*[0;1],xi,(i-1)/x_n,i/x_n);
    edge_dirichlet(i,1)=1;
    edge_dirichlet(i,2)=double(intg1);
     edge_dirichlet(edge_number-i+1,1)=1;
    edge_dirichlet(x_n*(y_n+1)+y_n*(x_n+1)-x_n+i,2)=double(intg2);
    
end
for i=1:y_n
    intg1=int(g4*[-1;0],eta,(i-1)/y_n,i/y_n);
    intg2=int(g3*[1;0],eta,(i-1)/y_n,i/y_n);
    edgedirichletnumber1=i*x_n+(i-1)*(x_n+1)+1;
    edge_dirichlet(edgedirichletnumber1,1)=1;
    edge_dirichlet(edgedirichletnumber1,2)=double(intg1);
     edge_dirichlet(edgedirichletnumber1+x_n,1)=1;
    edge_dirichlet(edgedirichletnumber1+x_n,2)=double(intg2);
end



element_boundary=zeros(element_number,1);
for i=1:x_n
    coordinate1=element_coordinate(i,1);
    x1=node_coordinate(coordinate1,1);
    y1=node_coordinate(coordinate1,2);
    element_boundary(i,1)=1;
    element_boundary(i,2)=function_p(x1+1/(2*x_n),y1+1/(2*y_n));

    coordinate1=element_coordinate(i+(y_n-1)*x_n,1);
    x1=node_coordinate(coordinate1,1);
    y1=node_coordinate(coordinate1,2);
    element_boundary(i+(y_n-1)*x_n,1)=1;
    element_boundary(i+(y_n-1)*x_n,2)=function_p(x1+1/(2*x_n),y1+1/(2*y_n));
end
for i=1:y_n
       coordinate1=element_coordinate((i-1)*x_n+1,1);
    x1=node_coordinate(coordinate1,1);
    y1=node_coordinate(coordinate1,2);
   element_boundary( (i-1)*x_n+1,1)=1;
 element_boundary( (i-1)*x_n+1,2)=function_p(x1+1/(2*x_n),y1+1/(2*y_n));

      coordinate1=element_coordinate(i*x_n,1);
    x1=node_coordinate(coordinate1,1);
    y1=node_coordinate(coordinate1,2);
    element_boundary( i*x_n,1)=1;
     element_boundary( i*x_n,1)=function_p(x1+1/(2*x_n),y1+1/(2*y_n));
end


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
% syms x y
% coefmat=coef();
% phi=[];
% phi1var=[ones(size(x)) ;x; y; x.*y ;y.*(1-x).*(1-y); x.*y.*(1-y)];
% phi2var=[ones(size(x)) ;x; y; x.*y ;x.*(1-x).*(1-y); x.*y.*(1-x)];
% for i=1:12
%     phi1=coefmat(i,1:6)*phi1var;
%     phi2=coefmat(i,7:12)*phi2var;
%     phicol=[phi1;phi2];
%     phi=[phi,phicol];
% end
% 
% phix=diff(phi,x);
% phiy=diff(phi,y);



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
load C11reference.mat;
load C12reference.mat;
load C22reference.mat;
load B11reference.mat;
load B12reference.mat;
load B22reference.mat;
load B21reference.mat;


load phi.mat;
load phix.mat;
load phiy.mat;

load phi1gaussvalue.mat;
load phi2gaussvalue.mat;

%%%%%%%%%%%%%%% gauss quadrature %%%%%%%%%%
gaussweight=[5/9,8/9,5/9];
gausspoint=[-sqrt(3/5),0,sqrt(3/5)];

xgausspoint2d=repmat(gausspoint,3,1);
ygausspoint2d=repmat(gausspoint.',1,3);
xgausspoint2d=reshape(xgausspoint2d,1,[]);
ygausspoint2d=reshape(ygausspoint2d,1,[]);

gaussweight2d=gaussweight.*gaussweight.';
gaussweight2d=reshape(gaussweight2d,1,[]);
gaussweight2d=repmat(gaussweight2d,localbasisfunctionnumber,1);

A=sparse(ubasis_function_number,ubasis_function_number);
B=sparse(ubasis_function_number,pbasis_function_number);
right=sparse(ubasis_function_number+pbasis_function_number,1);
imatrix=[];

%%%%%%%%%%%%gengerate  matrix A   B right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:element_number 
    local_coordinate1=node_coordinate(element_coordinate(i,1),1:2);
    local_coordinate2=node_coordinate(element_coordinate(i,2),1:2);
    local_coordinate3=node_coordinate(element_coordinate(i,3),1:2);
    local_coordinate4=node_coordinate(element_coordinate(i,4),1:2);
    J=[(local_coordinate2-local_coordinate1).'/2,(local_coordinate3-local_coordinate2).'/2];
    Fkb=(local_coordinate1+local_coordinate3).'/2;
    detJ=det(J);
    detJ=abs(detJ);
    invJ=inv(J);
    invJT=invJ.';
    C=invJ*invJT;
    xlocalgausspoint=J(1,1)*xgausspoint2d+J(1,2)*ygausspoint2d+Fkb(1,1)*ones(size(xgausspoint2d));
    ylocalgausspoint=J(2,1)*xgausspoint2d+J(2,2)*ygausspoint2d+Fkb(2,1)*ones(size(xgausspoint2d));
    
    %f1numerical=double(subs(f1,{xi,eta},{xlocalgausspoint,ylocalgausspoint}));
    %f1numerical=repmat(f1numerical,localbasisfunctionnumber,1);

     %f2numerical=double(subs(f2,{xi,eta},{xlocalgausspoint,ylocalgausspoint}));
    %f2numerical=repmat(f2numerical,localbasisfunctionnumber,1);
    fvalue=f(xlocalgausspoint,ylocalgausspoint);

    Alocal=nu*detJ*(C(1,1)*(A11reference+C11reference)+C(1,2)*(A12reference+A12reference.'+C12reference+C12reference.')+C(2,2)*(A22reference+C22reference));
    Blocal=detJ*(invJ(1,1)*B11reference+invJ(2,1)*B12reference+invJ(1,2)*B21reference+invJ(2,2)*B22reference);
     
    %rightlocal=detJ*(f1numerical.*phi1gaussvalue+f2numerical.*phi2gaussvalue).*gaussweight2d;
    rightlocal=detJ*(fvalue(1,:).*phi1gaussvalue+fvalue(2,:).*phi2gaussvalue).*gaussweight2d;
    rightlocal=sum(rightlocal,2);
    
    for j=1:localbasisfunctionnumber

        if j<=8
                jglobal2=mod(j-1,2)+1;
                jglobal1=element_coordinate(i,floor((j-1)/2)+1);
                jglobal=2*jglobal1+jglobal2-2;
               
            else
                jglobal1=j-8;
                jglobal=2*node_number+element_coordinate(i,jglobal1+4);
                
        end

        for k=1:localbasisfunctionnumber
            
            if k<=8
                kglobal2=mod(k-1,2)+1;
                kglobal1=element_coordinate(i,floor((k-1)/2)+1);
                kglobal=2*kglobal1+kglobal2-2;
               
            else
                kglobal1=k-8;
                kglobal=2*node_number+element_coordinate(i,kglobal1+4);
              
            end
            aa=[i;
            j;
            k;
            jglobal;
            kglobal];
            imatrix=[imatrix,aa];
            A(jglobal,kglobal)=A(jglobal,kglobal)+Alocal(j,k);
        end
        
            B(jglobal,i)=B(jglobal,i)+Blocal(j,1);
            right(jglobal,1)=right(jglobal,1)+rightlocal(j,1);
    end

  
end


Atotal=[A,B;B.',zeros(pbasis_function_number,pbasis_function_number)];
for i=1:node_number
    if(node_coordinate(i,3)==1)
        right=right-Atotal(:,2*i-1)*node_coordinate(i,4)-Atotal(:,2*i)*node_coordinate(i,5);
        Atotal(:,2*i)=0;
        Atotal(2*i,:)=0;
        Atotal(2*i,2*i)=1;
        Atotal(:,2*i-1)=0;
        Atotal(2*i-1,:)=0;
        Atotal(2*i-1,2*i-1)=1;
        right(2*i-1,1)=node_coordinate(i,4);
        right(2*i,1)=node_coordinate(i,5);
    end

end
% for i=1:node_number
%     if(node_coordinate(i,3)==1)
%         
% 
%         Atotal(2*i,:)=0;
%         Atotal(2*i,2*i)=1;
% 
%         Atotal(2*i-1,:)=0;
%         Atotal(2*i-1,2*i-1)=1;
%         right(2*i-1,1)=node_coordinate(i,4);
%         right(2*i,1)=node_coordinate(i,5);
%     end
% end
% 
for i=1:edge_number
    if(edge_dirichlet(i,1)==1)
        right=right-Atotal(:,2*node_number+i)*edge_dirichlet(i,2);
        Atotal(:,2*node_number+i)=0;
        Atotal(2*node_number+i,:)=0;
        Atotal(2*node_number+i,2*node_number+i)=1;
        right(2*node_number+i,1)=edge_dirichlet(i,2);
    end

end

 right=right-Atotal(:,ubasis_function_number+1)*element_boundary(1,2);
        Atotal(:,ubasis_function_number+1)=0;
        Atotal(ubasis_function_number+1,:)=0;
        Atotal(ubasis_function_number+1,ubasis_function_number+1)=1;
        right(ubasis_function_number+1,1)=element_boundary(1,2);
 
% for i=1:edge_number
%     if(edge_dirichlet(i,1)==1)
%         Atotal(2*node_number+i,:)=0;
%         Atotal(2*node_number+i,2*node_number+i)=1;
%         right(2*node_number+i,1)=edge_dirichlet(i,2);
%     end
% 
% end


basisfunctionweight=Atotal\right;


% xplot=0:(1/(3*x_n)):1;
% yplot=0:(1/(3*y_n)):1;
xplot=0:(1/(x_n)):1;
 yplot=0:(1/(y_n)):1;


[xplot,yplot]=meshgrid(xplot,yplot);
% xplot=[];
% yplot=[];
% for i=1:element_number
%      local_coordinate1=node_coordinate(element_coordinate(i,1),1:2);
%     local_coordinate2=node_coordinate(element_coordinate(i,2),1:2);
%     local_coordinate3=node_coordinate(element_coordinate(i,3),1:2);
%     local_coordinate4=node_coordinate(element_coordinate(i,4),1:2);
%     
% 
%     xplot=[xplot,local_coordinate1(1,1)+1/(3*x_n),local_coordinate1(1,1)+2/(3*x_n),local_coordinate1(1,1)+1/(3*x_n),local_coordinate1(1,1)+2/(3*x_n)];
%     yplot=[yplot,local_coordinate1(1,2)+1/(3*y_n),local_coordinate1(1,2)+1/(3*y_n),local_coordinate1(1,2)+2/(3*y_n),local_coordinate1(1,2)+2/(3*y_n)];
% end

% [uvalue,vvalue,pvalue]=plotvalue(xplot,yplot,x_n,y_n,ubasis_function_number,node_number,phi,element_coordinate,node_coordinate,basisfunctionweight);
uvalue=basisfunctionweight(1:2:node_number*2,1);
uvalue=reshape(uvalue,size(xplot));
vvalue=basisfunctionweight(2:2:node_number*2,1);
vvalue=reshape(vvalue,size(xplot));
pvalue=zeros(pbasis_function_number,1);
pvalue=basisfunctionweight(ubasis_function_number+1:ubasis_function_number+pbasis_function_number,1);
% pvalue=reshape(pvalue,size(xplot));
figure(1)
%quiver(xplot,yplot,uvalue,vvalue)
streamline(xplot,yplot,uvalue,vvalue)
%streamslice(xplot,yplot,uvalue,vvalue)
figure(2)
% plot3(xplot,yplot,pvalue)
colorbar
figure(3)
plot3(xplot,yplot,uvalue)
colorbar
figure(4)
plot3(xplot,yplot,vvalue)
colorbar
uvalueana=function_u(xplot,yplot);
vvalueana=function_v(xplot,yplot);
pvalueana=function_p(xplot,yplot);
figure(5)
plot3(xplot,yplot,pvalueana)
colorbar
figure(6)
plot3(xplot,yplot,uvalueana)
colorbar
figure(7)
plot3(xplot,yplot,vvalueana)
figure(8)
streamslice(xplot,yplot,uvalueana,vvalueana)
% % [xplotn,yplotn]=size(xplot);
% 
% syms x y
% phi11=matlabFunction(phi(1,1),'vars',[x y]);
% phi12=matlabFunction(phi(1,2),'vars',[x y]);
% phi13=matlabFunction(phi(1,3),'vars',[x y]);
% phi14=matlabFunction(phi(1,4),'vars',[x y]);
% phi15=matlabFunction(phi(1,5),'vars',[x y]);
% phi16=matlabFunction(phi(1,6),'vars',[x y]);
% phi17=matlabFunction(phi(1,7),'vars',[x y]);
% phi18=matlabFunction(phi(1,8),'vars',[x y]);
% phi19=matlabFunction(phi(1,9),'vars',[x y]);
% phi110=matlabFunction(phi(1,10),'vars',[x y]);
% phi111=matlabFunction(phi(1,11),'vars',[x y]);
% phi112=matlabFunction(phi(1,12),'vars',[x y]);
% phi21=matlabFunction(phi(2,1),'vars',[x y]);
% phi22=matlabFunction(phi(2,2),'vars',[x y]);
% phi23=matlabFunction(phi(2,3),'vars',[x y]);
% phi24=matlabFunction(phi(2,4),'vars',[x y]);
% phi25=matlabFunction(phi(2,5),'vars',[x y]);
% phi26=matlabFunction(phi(2,6),'vars',[x y]);
% phi27=matlabFunction(phi(2,7),'vars',[x y]);
% phi28=matlabFunction(phi(2,8),'vars',[x y]);
% phi29=matlabFunction(phi(2,9),'vars',[x y]);
% phi210=matlabFunction(phi(2,10),'vars',[x y]);
% phi211=matlabFunction(phi(2,11),'vars',[x y]);
% phi212=matlabFunction(phi(2,12),'vars',[x y]);
% for i=1:xplotn
%     for j=1:yplotn
%        xlocation=floor( xplot(i,j)*x_n);
%        ylocation=floor( yplot(i,j)*y_n);
%        if xlocation>=x_n
%            xlocation=x_n-1;
%        end
%        if ylocation>=y_n
%            ylocation=y_n-1;
%        end
%        elementlocation=xlocation+1+ylocation*x_n;
%        coordinatenumber1=element_coordinate(elementlocation,1);
%        coordinatenumber2=element_coordinate(elementlocation,2);
%        coordinatenumber3=element_coordinate(elementlocation,3);
%        coordinatenumber4=element_coordinate(elementlocation,4);
% 
%        local_coordinate1=node_coordinate(coordinatenumber1,1:2);
%     local_coordinate2=node_coordinate(coordinatenumber2,1:2);
%     local_coordinate3=node_coordinate(coordinatenumber3,1:2);
%     local_coordinate4=node_coordinate(coordinatenumber4,1:2);
%     edge1=element_coordinate(elementlocation,5);
%     edge2=element_coordinate(elementlocation,6);
%     edge3=element_coordinate(elementlocation,7);
%      edge4=element_coordinate(elementlocation,8);
%     J=[(local_coordinate2-local_coordinate1).'/2,(local_coordinate3-local_coordinate2).'/2];
%     Fkb=(local_coordinate1+local_coordinate3).'/2;
%     pvalue(i,j)=basisfunctionweight(ubasis_function_number+elementlocation,1);
%     coordinate=[xplot(i,j);yplot(i,j)];
%     local_coordinate=J*coordinate+Fkb;
%     weight=[basisfunctionweight(2*coordinatenumber1-1,1),basisfunctionweight(2*coordinatenumber1,1),...
%          basisfunctionweight(2*coordinatenumber2-1,1),basisfunctionweight(2*coordinatenumber2,1),...
%         basisfunctionweight(2*coordinatenumber3-1,1),basisfunctionweight(2*coordinatenumber3,1)...
%         basisfunctionweight(2*coordinatenumber4-1,1),basisfunctionweight(2*coordinatenumber4,1),...
%         basisfunctionweight(node_number*2+edge1,1),basisfunctionweight(node_number*2+edge2,1),...
%         basisfunctionweight(node_number*2+edge3,1),basisfunctionweight(node_number*2+edge4,1)];
%     phi1value=[phi11(coordinate(1,1),coordinate(2,1));phi12(coordinate(1,1),coordinate(2,1));
%         phi13(coordinate(1,1),coordinate(2,1));phi14(coordinate(1,1),coordinate(2,1));...
%         phi15(coordinate(1,1),coordinate(2,1));phi16(coordinate(1,1),coordinate(2,1));...
%         phi17(coordinate(1,1),coordinate(2,1));phi18(coordinate(1,1),coordinate(2,1));...
%         phi19(coordinate(1,1),coordinate(2,1));phi110(coordinate(1,1),coordinate(2,1));...
%         phi111(coordinate(1,1),coordinate(2,1));phi112(coordinate(1,1),coordinate(2,1))];
%      phi2value=[phi21(coordinate(1,1),coordinate(2,1));phi22(coordinate(1,1),coordinate(2,1));...
%         phi23(coordinate(1,1),coordinate(2,1));phi24(coordinate(1,1),coordinate(2,1));...
%         phi25(coordinate(1,1),coordinate(2,1));phi26(coordinate(1,1),coordinate(2,1));...
%         phi27(coordinate(1,1),coordinate(2,1));phi28(coordinate(1,1),coordinate(2,1));...
%         phi29(coordinate(1,1),coordinate(2,1));phi210(coordinate(1,1),coordinate(2,1));...
%         phi211(coordinate(1,1),coordinate(2,1));phi212(coordinate(1,1),coordinate(2,1))];
%      uvalue(i,j)=weight*phi1value;
%      vvalue(i,j)=weight*phi2value;
%       
%     end
% end
