function [msk]=genMsk2(GVM)
% remove noisy data
[s1,s2,~]=size(GVM);
tmp=(GVM(:,:,8)>2.5);
CNR=reshape(tmp,[s1 s2]);
tmp=(GVM(:,:,9)>2.5);
CNT=reshape(tmp,[s1 s2]);
tmp=(GVM(:,:,6)>0);
MSK=reshape(tmp,[s1 s2]);
msk=MSK & CNR & CNT;