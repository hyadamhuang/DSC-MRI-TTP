function [MSK]=findVOLMSK(VOL,flag)
% input: VOL(x,y,t,z), 3D MRP data, t temporal samples
% output: MSK(x,y,1,z), brain tissue mask
[s1,s2,~,s4]=size(VOL);
MSK=(VOL(:,:,1,:)>0);
se=strel('disk',1);
for z=1:s4
    I=double(VOL(:,:,1,z));
    mx=max(I(:));
    mn=min(I(:));
    level = graythresh(I);
    BW = im2bw((I-mn)/(mx-mn),level/5);
    BW=imopen(BW,se);
    CC=bwconncomp(BW);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    for i=1:length(numPixels)
        if numPixels(i)<100
            BW(CC.PixelIdxList{i})=0;
        end
    end
    BW=imclose(BW,se);
    MSK(:,:,1,z)=BW;
end
% find the largest connected component as the brain tissue mask
CC=bwconncomp(MSK);%reshape(MSK,s1,s2,s4));
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx]=max(numPixels);
MSK(:)=0;
MSK(CC.PixelIdxList{idx})=1;
if flag>0
    montage(reshape(MSK,s1,s2,1,s4),'Size',[ceil(s4/6) NaN]);
end
    