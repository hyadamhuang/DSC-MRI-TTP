function [GVM,MSK]=compGVMii(VOL1,MSK,dt)
% input: VOL1(x,y,t), 2D MRP data, MSK(x,y) mask, dt sampling time interval
% output: GVM(x,y,10) model, MSK(x,y), updated brain tissue mask
[s1,s2,s3]=size(VOL1);%x,y,t
GVM=zeros(s1,s2,10);%10 parameters: BAT,TTP,T1,T2,GMX,alpha,AUC,CNRi,CNRx,ds(=abs(S(t=BAT)-S(t=TTP)))
% GVM's 10 parameters: refer to: Huang A, Lee CW, Liu HM. Curve fitting
% criteria to determine arterial input function for MR perfusion analysis.
% Venice, Italy: Proceedings of the 16th IEEE International Symposium on
% Biomedical Imaging; 2019

for i=1:s1
    for j=1:s2
        if MSK(i,j)>0
            S=double(squeeze(VOL1(i,j,:)));
            [gvm,msk]=compOneVoxel(S,dt);
            MSK(i,j)=msk;
            for k=1:10
                GVM(i,j,k)=gvm(k);
            end
        end
    end
end
