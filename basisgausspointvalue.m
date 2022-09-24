function [phi1gaussvalue,phi2gaussvalue]=basisgausspointvalue(phi)

syms x y;
gaussweight=[5/9,8/9,5/9];
gausspoint=[-sqrt(3/5),0,sqrt(3/5)];

xgausspoint2d=repmat(gausspoint,3,1);
ygausspoint2d=repmat(gausspoint.',1,3);
xgausspoint2d=reshape(xgausspoint2d,1,[]);
ygausspoint2d=reshape(ygausspoint2d,1,[]);

 phi2gaussvalue=[];
  phi1gaussvalue=[];
for i=1:12
    phi1=double(subs(phi(1,i),{x,y},{xgausspoint2d,ygausspoint2d}));
    phi2=double(subs(phi(2,i),{x,y},{xgausspoint2d,ygausspoint2d}));

    phi1gaussvalue=[phi1gaussvalue;phi1];
    phi2gaussvalue=[phi2gaussvalue;phi2];
end

save phi1gaussvalue.mat phi1gaussvalue;
save phi2gaussvalue.mat phi2gaussvalue;