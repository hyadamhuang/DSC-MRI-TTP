% % % % % %
% This code is a part of manuscript
% Run this program to draw Figure 1
% This code has been submitted as supplementary for publication review
% Journal: Scientific Reports
% Manuscript title: "Time to Peak and Full Width at Half Maximum in MR Perfusion: Valuable Indicators for Monitoring Moyamoya Patients after Revascularization"
% Authors: Adam Huang (National Central University, Taiwan), Chung-Wei Lee
% (National Taiwan University Hospital), Hon-Man Liu (Fu Jen Catholic
% University Hospital, Fu-Jen Catholic University, New Taipei City, Taiwan)
% Date: 2020/11/20
% % % % % %
%% read data
% read DSC-MRI data at convexity VOL15, 128x128x50, 50 temporal samples, dt=2s
load('convexity.mat','VOL15');
% read DSC-MRI data at cerebellum VOL6, 128x128x50, 50 temporal samples, dt=2s
load('cerebellum.mat','VOL6');
ref=[84 55;82 78]; % cerebellum references (right & left)
thres=[2,4,6]; % thresholds
x=66;y=96; % randomly chosen sample point at convexity location (66,96)
%% GVM computation for deriving TTP and FWHM map
[TTP3,FWHM,BAR,tA,fA]=compTtpFwhm2D(VOL15,VOL6,ref,thres);
BAR=imresize(BAR,[512 50]); % resize color bar
TTP3=imresize(TTP3,[512 512]); % resize TTP map
FWHM=imresize(FWHM,[512 512]); % resize FWHMn map
%% draw to Figure 1 (lower-left, BAR, TTP, FWHMn)
% create a large area to draw Figure 1
I=ones(3*600,6*600,3,'uint8')*255;
I((1:512)+1180-50,1:50,:)=BAR;
I((1:512)+1180-50,(1:512)+200,:)=TTP3;
I((1:512)+1180-50,(1:512)+850,:)=FWHM;
%% draw Figure 1 (upper-left, DSC-MRI images at convexity  and cerebellum)
J=zeros(128,128,3,'uint8');
n=[11 13 15 17]; % 4 sampling time points
% convexity
for i=1:4
    II=double(VOL15(:,:,n(i),1));% MRI image at convexity, time n(i)
    mx=max(II(:));
    mn=min(II(:));
    II=uint8((II-mn)/(mx-mn)*255);% adjust image brightness
    J(:,:,1)=II; % red channel
    J(:,:,2)=II; % green channel
    J(:,:,3)=II; % blue channel
    % mark sample point at x=66;y=96; in a green box
    for ii=x-2:x+2
        for jj=y-2:y+2
            if (abs(ii-x)>1) || (abs(jj-y)>1)
                J(ii,jj,1)=0;
                J(ii,jj,2)=255;
                J(ii,jj,3)=0;
            end
        end
    end
    % resize image
    JJ=imresize(J,[384,384]);
    % draw to I (Figure 1)
    I((1:384)+80,(1:384)+(i-1)*400,:)=JJ;
end
% cerebellum
for i=1:4
    II=double(VOL6(:,:,n(i),1)); % MRI image at cerebellum, time n(i)
    mx=max(II(:));
    mn=min(II(:));
    II=uint8((II-mn)/(mx-mn)*255);
    J(:,:,1)=II; % red channel
    J(:,:,2)=II; % green channel
    J(:,:,3)=II; % blue channel
    % mark reference points at ref in green boxes
    for kk=1:size(ref,1)
        for ii=ref(kk,1)-5:ref(kk,1)+5
            for jj=ref(kk,2)-5:ref(kk,2)+5
                if (abs(ii-ref(kk,1))>4) || (abs(jj-ref(kk,2))>4)
                    J(ii,jj,1)=0;
                    J(ii,jj,2)=255;
                    J(ii,jj,3)=0;
                end
            end
        end
    end
    JJ=imresize(J,[384,384]);
    I((1:384)+540,(1:384)+(i-1)*400,:)=JJ;
end
%% draw figure 1 (image part)
figure,imshow(I);
hold on,
%% draw Figure (Right-side, GVM computation)
% raw DSC-MRI signal S(t)
S=double(squeeze(VOL15(x,y,:,1))); % example sample point at (x,y)
ix=(0:length(S)-1)/50*1700+1850;
t=0:2:(length(S)-1)*2;% scanning time interval = 2 seconds
S(1)=mean(S(2:5));% correcting the first image (for NTUH scanning protocol)
S1=max(S);
S0=min(S);
SS=(S-S0)/(S1-S0);
plot(ix,-SS*400+600,'k.'); % draw S(t)
plot(ix(n),-SS(n)*400+600,'gs');% 4 time points of the left MRI images
text(1860,300,'S(t)');
% signal to time concentration curve, LinearLinearModel to find BAT
C=-log(S/max(S));
[BAT,TTP,Cpre,CC,CMX,CNR,YY]=BAT_LLM(t',C,0);% fit LinearLinearModel
plot(ix,-C*100+1000,'k.');% C(t)
plot(ix(1:length(YY)),-C(1:length(YY))*100+1000,'ro');% fit samples used
plot(ix(1:length(YY)),-YY*100+1000,'g');% LinearLinearModel
text(1860,1100,'C(t)');
tmp=diff(YY);
i=find(tmp>0);
text(ix(i(1))-15,-YY(i(1))*100+1000+100-25,'t0');
% G(t): Rigid GVM
[TTP,k1,k2,alpha,beta,K,YMX]=findTTP_LGV_i(t',C,BAT,TTP,0.5,0);
if alpha>0
    tt=t(1):0.1:t(end);% finer time interval
    ttp=TTP-BAT;
    u=tt-BAT;
    IX=(u<=0);
    u(IX)=0;
    y1=K*u.^alpha.*exp(-u/beta);%exp(B1)/ttp^alpha*exp(alpha)*(u).^alpha.*exp(-alpha/ttp*(u));% from data model
end
plot(tt/100*1700+1850,-y1*100+1500,'m',ix,-C*100+1500,'k.',ix(k1:k2),-C(k1:k2)*100+1500,'ro');
[mx,i]=max(y1);
text(tt(i)/100*1700+1850,-y1(i)*100+1500-50,'TTP0');
i=find(y1>0.5*mx);
text(tt(i(end))/100*1700+1850+20,-y1(i(end))*100+1500-15,'FWHM = t2-t1');
y_h=(y1(i(1))+y1(i(end)))/2;
plot([tt(i(1)) tt(i(1)) tt(i(end)) tt(i(end))]/100*1700+1850,...
    [0 -y_h -y_h 0]*100+1500,'b');
text(tt(i(1))/100*1700+1850-15,1500+15,'t1');
text(tt(i(end))/100*1700+1850-15,1500+15,'t2');
text(1860,1580,'G(t)');
text(2060,1610,'TTP=TTP0-TTPc+1');
text(2060,1700,'FWHMn=10(FWHM-FWHMc)/FWHMc+1');
%% other text notes
plot([1801 3600 3600 1801 1801],[1 1 1800 1800 1],'k','LineWidth',1);
text(1900,50,'Rigid GVM Computation');

text(50,40,'Convexity (DSC-MRI)');
text(50,500,'Cerebellum');
text(1600,270,'* * *');
text(1600,750,'* * *');

text(250,1080,'TTP');
text(900,1080,'FWHMn');
text(60,1680-50,'0s');
text(60,1430-50,'6s');
text(60,1180-50,'12s');

s1=[];
for i=1:length(tA)-1
    s1=[s1,num2str(round(tA(i))),'/'];
end
s1=[s1,num2str(round(tA(end))),'% ('];
for i=1:length(thres)-1
    s1=[s1,num2str(round(thres(i))),'/'];
end
s1=[s1,num2str(round(thres(end))),'s)'];
text(200,1700,s1);
s2=[];
for i=1:length(fA)-1
    s2=[s2,num2str(round(fA(i))),'/'];
end
s2=[s2,num2str(round(fA(end))),'% ('];
for i=1:length(thres)-1
    s2=[s2,num2str(round(thres(i))),'/'];
end
s2=[s2,num2str(round(thres(end))),'s)'];
text(850,1700,s2);
%% end of draw
hold off;
