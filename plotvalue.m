function [uvalue,vvalue,pvalue]=plotvalue(xplot,yplot,x_n,y_n,ubasis_function_number,node_number,phi,element_coordinate,node_coordinate,basisfunctionweight)

[xplotn,yplotn]=size(xplot);

syms x y
phi11=matlabFunction(phi(1,1),'vars',[x y]);
phi12=matlabFunction(phi(1,2),'vars',[x y]);
phi13=matlabFunction(phi(1,3),'vars',[x y]);
phi14=matlabFunction(phi(1,4),'vars',[x y]);
phi15=matlabFunction(phi(1,5),'vars',[x y]);
phi16=matlabFunction(phi(1,6),'vars',[x y]);
phi17=matlabFunction(phi(1,7),'vars',[x y]);
phi18=matlabFunction(phi(1,8),'vars',[x y]);
phi19=matlabFunction(phi(1,9),'vars',[x y]);
phi110=matlabFunction(phi(1,10),'vars',[x y]);
phi111=matlabFunction(phi(1,11),'vars',[x y]);
phi112=matlabFunction(phi(1,12),'vars',[x y]);
phi21=matlabFunction(phi(2,1),'vars',[x y]);
phi22=matlabFunction(phi(2,2),'vars',[x y]);
phi23=matlabFunction(phi(2,3),'vars',[x y]);
phi24=matlabFunction(phi(2,4),'vars',[x y]);
phi25=matlabFunction(phi(2,5),'vars',[x y]);
phi26=matlabFunction(phi(2,6),'vars',[x y]);
phi27=matlabFunction(phi(2,7),'vars',[x y]);
phi28=matlabFunction(phi(2,8),'vars',[x y]);
phi29=matlabFunction(phi(2,9),'vars',[x y]);
phi210=matlabFunction(phi(2,10),'vars',[x y]);
phi211=matlabFunction(phi(2,11),'vars',[x y]);
phi212=matlabFunction(phi(2,12),'vars',[x y]);
for i=1:xplotn
    for j=1:yplotn
       xlocation=floor( xplot(i,j)*x_n);
       ylocation=floor( yplot(i,j)*y_n);
       if xlocation>=x_n
           xlocation=x_n-1;
       end
       if ylocation>=y_n
           ylocation=y_n-1;
       end
       elementlocation=xlocation+1+ylocation*x_n;
       coordinatenumber1=element_coordinate(elementlocation,1);
       coordinatenumber2=element_coordinate(elementlocation,2);
       coordinatenumber3=element_coordinate(elementlocation,3);
       coordinatenumber4=element_coordinate(elementlocation,4);

       local_coordinate1=node_coordinate(coordinatenumber1,1:2);
    local_coordinate2=node_coordinate(coordinatenumber2,1:2);
    local_coordinate3=node_coordinate(coordinatenumber3,1:2);
    local_coordinate4=node_coordinate(coordinatenumber4,1:2);
    edge1=element_coordinate(elementlocation,5);
    edge2=element_coordinate(elementlocation,6);
    edge3=element_coordinate(elementlocation,7);
     edge4=element_coordinate(elementlocation,8);
    J=[(local_coordinate2-local_coordinate1).'/2,(local_coordinate3-local_coordinate2).'/2];
    Fkb=(local_coordinate1+local_coordinate3).'/2;
    pvalue(i,j)=basisfunctionweight(ubasis_function_number+elementlocation,1);
    coordinate=[xplot(i,j);yplot(i,j)];
    local_coordinate=J*coordinate+Fkb;
    weight=[basisfunctionweight(2*coordinatenumber1-1,1),basisfunctionweight(2*coordinatenumber1,1),...
         basisfunctionweight(2*coordinatenumber2-1,1),basisfunctionweight(2*coordinatenumber2,1),...
        basisfunctionweight(2*coordinatenumber3-1,1),basisfunctionweight(2*coordinatenumber3,1)...
        basisfunctionweight(2*coordinatenumber4-1,1),basisfunctionweight(2*coordinatenumber4,1),...
        basisfunctionweight(node_number*2+edge1,1),basisfunctionweight(node_number*2+edge2,1),...
        basisfunctionweight(node_number*2+edge3,1),basisfunctionweight(node_number*2+edge4,1)];
    phi1value=[phi11(coordinate(1,1),coordinate(2,1));phi12(coordinate(1,1),coordinate(2,1));
        phi13(coordinate(1,1),coordinate(2,1));phi14(coordinate(1,1),coordinate(2,1));...
        phi15(coordinate(1,1),coordinate(2,1));phi16(coordinate(1,1),coordinate(2,1));...
        phi17(coordinate(1,1),coordinate(2,1));phi18(coordinate(1,1),coordinate(2,1));...
        phi19(coordinate(1,1),coordinate(2,1));phi110(coordinate(1,1),coordinate(2,1));...
        phi111(coordinate(1,1),coordinate(2,1));phi112(coordinate(1,1),coordinate(2,1))];
     phi2value=[phi21(coordinate(1,1),coordinate(2,1));phi22(coordinate(1,1),coordinate(2,1));...
        phi23(coordinate(1,1),coordinate(2,1));phi24(coordinate(1,1),coordinate(2,1));...
        phi25(coordinate(1,1),coordinate(2,1));phi26(coordinate(1,1),coordinate(2,1));...
        phi27(coordinate(1,1),coordinate(2,1));phi28(coordinate(1,1),coordinate(2,1));...
        phi29(coordinate(1,1),coordinate(2,1));phi210(coordinate(1,1),coordinate(2,1));...
        phi211(coordinate(1,1),coordinate(2,1));phi212(coordinate(1,1),coordinate(2,1))];
     uvalue(i,j)=weight*phi1value;
     vvalue(i,j)=weight*phi2value;
      
    end
end

end