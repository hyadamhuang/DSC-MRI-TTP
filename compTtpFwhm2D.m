function [TTP,FWHM,BAR,tA,fA]=compTtpFwhm2D(VOL15,VOL6,ref,thres)
% input:
%   VOL15(x,y,t), 2D DSC-MRI data at convexity, t temporal samples, assuming dt(sampling time)=2s
%   VOL6(x,y,t) 2D DSC-MRI data at cerebellum, t temporal samples, assuming dt(sampling time)=2s
%   ref, coordinates of references at cerebellum
%   thres, time thresholds (2s, 4s, 6s)
% output:
%   TTP, TTP color map
%   FWHM, FWHM color map
%   BAR, BAR color map (0~12s)
%   tA, % of area at the given thres for TTP
%   fA, % of area at the given thres for FWHM
[s1,s2,s3]=size(VOL15);
[MSK15]=findVOLMSK(VOL15,0);
[GVM15,MSK15]=compGVMii(VOL15,MSK15,2);
[MSK6]=findVOLMSK(VOL6,0);
[GVM6,MSK6]=compGVMii(VOL6,MSK6,2);
len=size(ref,1);
if len<1 % need at least one reference
    TTP=[];FWHM=[];BAR=[];tA=[];fA=[];
    return;
end
ttp0=0;
ttp6=0;
for k=1:size(ref,1) 
    aa=ref(k,1);bb=ref(k,2);
    ttp2=0;nt=0;ttp4=0;
    for a=aa-4:aa+4
        for b=bb-4:bb+4
            if MSK6(a,b)>0
                ttp2=ttp2+GVM6(a,b,2,1); %ttp
                ttp4=ttp4+(GVM6(a,b,4,1)-GVM6(a,b,3,1));%t2-t1
                nt=nt+1;
            end
        end
    end
    ttp0=ttp0+ttp2/nt;
    ttp6=ttp6+ttp4/nt;
end
ttp0=ttp0/len;% reference TTPc
ttp6=ttp6/len;% reference FWHMc
% time map color
color=zeros(24,3);
color(1,:)=[100,51,255];
color(2,:)=[75,51,255];
color(3,:)=[51,51,255];
color(4,:)=[51,100,255];
color(5,:)=[51,150,255];
color(6,:)=[51,200,255];
color(7,:)=[51,255,255];
color(8,:)=[51,255,175];
color(9,:)=[51,255,91];
color(10,:)=[122,255,71];
color(11,:)=[193,255,51];% 6
color(12,:)=[213,245,54];
color(13,:)=[236,236,57];% 7
color(14,:)=[245,215,30];
color(15,:)=[255,195,0];% 8
color(16,:)=[255,120,25];
color(17,:)=[255,87,51];% 9
color(18,:)=[255,45,25];
color(19,:)=[255,0,0];% 10
color(20,:)=[200,0,0];
color(21,:)=[150,0,0];% 11
color(22,:)=[125,0,0];
color(23,:)=[100,0,0];% 12
color(24,:)=[75,0,0];% 12

[msk15]=genMsk2(GVM15);

    
BAR=zeros(240,10,3,'uint8');% color bar
for abc=1:24
    BAR((abc-1)*10+1:abc*10,:,1)=color(25-abc,1);
    BAR((abc-1)*10+1:abc*10,:,2)=color(25-abc,2);
    BAR((abc-1)*10+1:abc*10,:,3)=color(25-abc,3);
end
tnum=length(thres);
tA=zeros(tnum,1);
fA=zeros(tnum,1);
TTP=zeros(128,128,3,'uint8');
for i=1:s1
    for j=1:s2
        if msk15(i,j)>0
            t=GVM15(i,j,2);
            ttpn=t-ttp0+1;
            for k=1:tnum
                if ttpn<thres(k)
                    tA(k)=tA(k)+1;
                end
            end
            jx=min(24,max(1,floor(ttpn/0.5)+1));% 0.5 sec per color
            TTP(i,j,:)=color(jx,:);
        end
    end
end
FWHM=zeros(128,128,3,'uint8');
for i=1:s1
    for j=1:s2
        if msk15(i,j)>0
            t=GVM15(i,j,4)-GVM15(i,j,3);
            fwhmn=10*(t-ttp6)/ttp6+1;
            for k=1:tnum
                if fwhmn<thres(k)
                    fA(k)=fA(k)+1;
                end
            end
            jx=min(24,max(1,floor(fwhmn/0.5)+1));% 0.5 sec per color
            FWHM(i,j,:)=color(jx,:);
        end
    end
end
T=(msk15>0);
tot=sum(T(:));
fA=fA/tot*100;
tA=tA/tot*100;