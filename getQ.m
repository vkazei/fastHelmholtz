function Q = getQ(h,n,zt,xt)
% same as getP but with normalization
% that mimics the delta function
Q = getP(h,n,zt,xt)/prod(h);
