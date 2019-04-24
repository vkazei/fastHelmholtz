function Q = getQ(zt,xt,z,x)
% same as getP but with normalization
% that mimics the delta function
Q = getP(zt,xt,z,x)*(x(1)-x(2))^-2;
