function [phi,phix,phiy]=referencebasicfunction(coefmat)

syms x y
phi=[];
phi1var=[1 ;x; y; x.*y ;-(y+1).*(y-1).*(x-1);(x+1).*(y-1).*(y+1)];
phi2var=[1 ;x; y; x.*y ;-(x+1).*(x-1).*(y-1);(x+1).*(x-1).*(y+1)];
for i=1:12
    phi1=coefmat(i,1:6)*phi1var;
    phi2=coefmat(i,7:12)*phi2var;
    phicol=[phi1;phi2];
    phi=[phi,phicol];
end

phix=diff(phi,x);
phiy=diff(phi,y);

save phix.mat phix;
save phiy.mat phiy;
save phi.mat phi;
end