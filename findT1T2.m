function [T1,T2]=findT1T2(K,alpha,beta)
u=0:0.1:60;
c=K*u.^alpha.*exp(-u/beta);
ix=(c>=(0.5*max(c)));
T2=max(u(ix));
T1=min(u(ix));