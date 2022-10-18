function fvalue=f(x,y)

% fvalue(1,:)=0*ones(size(x));%sin(x);%ones(size(x));
% fvalue(2,:)=0*ones(size(x));%sin(y);%ones(size(y));
%fvalue(1,:)=(2*pi*cos(2*pi*x)+pi^2*sin(pi*x)+pi^2*sin(pi*y));
%fvalue(2,:)=(2*pi*cos(2*pi*y)-pi^3*sin(pi*x).*y);
nu=1;
%  fvalue(1,:)=-2*nu*(x.^2+y.^2)-nu*exp(-y)+pi^2*cos(pi*x).*cos(2*pi*y);
% fvalue(2,:)=- 4*nu*x.*y-nu.*pi^3*sin(pi*x)+2*pi*(2-pi*sin(pi*x)).*sin(2*pi*y);
 fvalue(1,:)=(1-4*nu)*ones(size(x));
fvalue(2,:)=(1-4*nu)*ones(size(x));
end