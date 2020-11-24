function [BAT,TTP,Cpre,CC,CMX,CNR,YY]=BAT_LLM(t,y,flag)
% original paper: An automatic approach for estimating bolus arrival time
% in dynamic contrast MRI using piecewise continuous regression models
% by LH Cheong, TS Koh, and Z Hou
% input:
%   t: time
%   y: contrast intensity
%   flag: drawing flag
% ouput:
%   BAT: balus arrival time
%   TTP: time to peak
%   Cpre: baseline
%   CC: y-Cpre
%   CMX: max contrast
%   CNR: contrast intensity noise ratio (t:0~BAT)
% Author: Adam Huang, 2020/07/24
% adamhuan@gmail.com
[~,len]=max(y); % find the peak time position
len=max(5,len); % TTP occurs at least >= t(5) 
TTP=t(len);
[c1,c2]=size(y);
if c2>c1
    y=y';
end
C=y(1:len);
ix=1;
err=realmax;% a big number
if flag>0
    figure,
    plot(t,y,'r');
end
for i=2:len-1
    X=ones(len,2);
    X(:,2)=0;
    for j=i+1:len
        X(j,2)=t(j)-t(i);
    end
%     B=inv(transpose(X)*X)*(transpose(X)*C); % find inv() explicitly
    B=(transpose(X)*X)\(transpose(X)*C); % '\' using gauss elimination
    D=C-X*B;
    e=norm(D);
    if e<err
        ix=i;
        err=e;
        if flag>0
            hold on,
            plot(t(1:len),X*B,'g');
            hold off;
        end
    end
end
BAT=t(ix);
Cpre=mean(y(1:ix));
CC=y-Cpre;
CMX=max(CC(:));
sd=std(CC(1:ix)-Cpre);
if sd<=eps
    CNR=0;
else
    CNR=CMX/sd;
end
% % draw result if flag is set
X=ones(len,2);
X(:,2)=0;
for j=ix:len
    X(j,2)=t(j)-BAT;
end
%     B=inv(transpose(X)*X)*(transpose(X)*C);
B=(transpose(X)*X)\(transpose(X)*C);
YY=X*B;
if flag>0
    figure,plot(t,CC,'r.',t(1:len),YY-Cpre,'b');
end