function fvalue=f(x,y)

% fvalue(1,:)=0*ones(size(x));%sin(x);%ones(size(x));
% fvalue(2,:)=0*ones(size(x));%sin(y);%ones(size(y));
fvalue(1,:)=(2*pi*cos(2*pi*x)+pi^2*sin(pi*x)+pi^2*sin(pi*y));
fvalue(2,:)=(2*pi*cos(2*pi*y)-pi^3*sin(pi*x).*y);
end