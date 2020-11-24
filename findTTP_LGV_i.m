function [TTP,k1,k2,alpha,beta,K,YMX]=findTTP_LGV_i(t,C,BAT,TTP,q2,flag)
% k1 the first non-zero gamma variate function point
% k2 defines the 50% drop point of gamma variate function
% fit Madsen's linear model to C(k1:k2) 
% % step 0: preprocessing
[c1,c2]=size(t);
if c2>c1
    t=t';
end
[c1,c2]=size(C);
if c2>c1
    C=C';
end
% % step 1: find k1 (gamma function nonzero starting point) 
dt=t(2)-t(1);
k1=floor(BAT/dt)+2;
% % step 2: find k2 (gamma function drop to half)
len=length(C);
CMX=max(C);
c2=CMX*q2;
ix=min(len-5,floor(TTP/dt)+1);
for i=ix+2:len
    if C(i)<=c2
        k2=i;
        break;
    end
    k2=i;
end
if k2==len % very long duration
    k2=floor((len-ix)/2)+ix;
elseif k2-k1<3 % very short duration
    k2=k2+1;
end
% % step 3: find ALPHA,CMX,TTP by Madsen's linear fitting to Gamma Variate
x=t(k1:k2);
y=C(k1:k2);
if flag>0
    figure,
    plot(t,C,'r.',x,y,'bo');
end
ei=realmax;
alpha=0;beta=0;K=0;YMX=0;
for i=1:2
    if i<2
        x0=TTP-(TTP-BAT)/2;
        x1=min(t(len-5),TTP+6);
        dx=0.5; % linear search with a time step of 0.5 sec
    else
        x0=TTP-min(0.4,(TTP-BAT)/2);
        x1=TTP+0.4;
        dx=0.1; % linear search with a time step of 0.1 sec
    end
    for ttpx=x0:dx:x1 % linear search
        ttp=ttpx-BAT;
        if ttp>0
            xxt=(x-BAT)/ttp;
            IX=((xxt>0.01)&(y>0.01*CMX));
            N=sum(IX);
            if N>=3
                X=ones(N,2);
                X(:,2)=1+log(xxt(IX))-xxt(IX);
                Ci=log(y(IX));
%                 B=inv(transpose(X)*X)*(transpose(X)*Ci);
                B=(transpose(X)*X)\(transpose(X)*Ci);
                a=B(2);
                u=t(k1:k2)-BAT;
                IX=(u<=0);
                u(IX)=0;
                y1=exp(B(1))/ttp^a*exp(a)*(u).^a.*exp(-a/ttp*(u));% from data model
                e=norm(y1-y);
                if e<ei
                    ei=e;
                    TTP=ttpx;
                    B1=B(1);
                    alpha=a;
                    if flag>0
                        hold on,
                        plot(u+BAT,y1,'g.');
                        hold off;
                    end
                end
            end
        end
    end
end
if alpha>0
    ttp=TTP-BAT;
    beta=ttp/alpha;
    K=exp(B1)/ttp^alpha*exp(alpha);
    YMX=exp(B1);
    if flag>0
        u=t-BAT;
        IX=(u<=0);
        u(IX)=0;
        y1=K*u.^alpha.*exp(-u/beta);%exp(B1)/ttp^alpha*exp(alpha)*(u).^alpha.*exp(-alpha/ttp*(u));% from data model
        hold on,plot(t,y1,'m',t,C-y1,'g');hold off;
    end
end