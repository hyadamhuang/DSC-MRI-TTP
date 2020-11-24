function [GVM,MSK]=compOneVoxel(S,dt)
% rigid gamma variate modelling
% input:
%   S: signal intensity
%   dt: scanning time interval
% ouput:
%   GVM: model parameters
%   MSK: brain tissue mask
% Author: Adam Huang, 2020/07/24
% adamhuan@gmail.com
len=length(S);
T=0:dt:dt*(len-1);
T=T';
S=double(S);
mn=min(S);
GVM=zeros(10,1);
MSK=1;
% % parameter 1 % %
if mn<=50 % signal near border with too small intensity
    GVM(:)=0;
else % mn>50 % signal value cannot be too small
    S(1)=mean(S(2:5)); % first point protocol correction for NTUH
    C=-log(S/max(S)); % signal -> contrast concentration
    [BAT,TTP,Cpre,C,CMX,CNR,~]=BAT_LLM(T,C,0);
    % % parameter 2 % %
    if CNR<4 % CMX/std: too noisy or CMX too small
        MSK(:)=0;
    else % if CNR>=4
        [TTP,k1,k2,alpha,beta,K,GMX]=findTTP_LGV_i(T,C,BAT,TTP,0.5,0);
        % % parameter 3, fit successfully % %
        if alpha>0
            [T1,T2]=findT1T2(K,alpha,beta); % mid-height
            k2=min(ceil(T2/dt)+1,len);
            u=T-BAT;
            IX=(u<=0);
            u(IX)=0;
            Y=K*u.^alpha.*exp(-u/beta);
            CNRi=GMX/std(C(1:k2)-Y(1:k2));
            CNRx=GMX/std(C-Y);
            AUC=trapz(T,[C(1:k2-1);Y(k2:end)]);
%             try
%                 AUC=trapz(T(1:k2),C(1:k2))+trapz(T(k2:end),Y(k2:end));
%             catch
%                 disp(['k2=',num2str(k2)]);
%             end
%             if CNR>=4
            GVM(1)=BAT;
            GVM(2)=TTP;
            GVM(3)=T1;
            GVM(4)=T2;
            GVM(5)=GMX;
            GVM(6)=alpha;
            GVM(7)=AUC;
            GVM(8)=CNRi;
            GVM(9)=CNRx;
            % abs(S(t=BAT)-S(t=TTP))
            GVM(10)=abs(S(floor(BAT/dt))-min(S(floor(TTP/dt)),S(floor(TTP/dt)+1)));
%             end
        end
    end
end