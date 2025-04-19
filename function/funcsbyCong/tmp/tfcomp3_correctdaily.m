function [corrspec1_2,corrspec1_32,corrspec1_432,corrtime1_2,corrtime1_32,corrtime1_432,...
    corrtime1_2_day,corrtime1_32_day,corrtime1_432_day,corrgood1_2,corrgood1_32,corrgood1_432,label_list,...
    noise1_2,noise3_2,noise1_32,noise4_32,noise1_432] =...
    tfcomp3_correctdaily(TF_cal,TFs,spec_mat,Zraw,H1raw,H2raw,Praw,Zsig,H1sig,H2sig,Psig,NFFT,dt,f,npad0,npts,goodP,goodH,goodZ,taxisZ,taptime)

TF_name = TF_cal{1};

comp1 = TF_name(1);
comp2 = TF_name(5);
comp3 = TF_name(4);
comp4 = TF_name(2);

TFidx12 = strmatch([comp1,comp2],{TFs.label},'exact');
TFidx32 = strmatch([comp3,comp2],{TFs.label},'exact');
TFidx42 = strmatch([comp4,comp2],{TFs.label},'exact');
TFidx43_2 = strmatch([comp4,comp3,'-',comp2],{TFs.label},'exact');
TFidx13_2 = strmatch([comp1,comp3,'-',comp2],{TFs.label},'exact');
TFidx14_32 = strmatch([comp1,comp4,'-',comp3,comp2],{TFs.label},'exact');

label_list{1} = [comp1,comp2];
label_list{2} = [comp1,comp3,'-',comp2];
label_list{3} = [comp1,comp4,'-',comp3,comp2];

if strcmp(comp1,'Z')==1
    spec_1 = spec_mat(:,1);
    isgood1 = goodZ;
elseif strcmp(comp1,'1')==1
    spec_1 = spec_mat(:,2);
    isgood1 = goodH;
elseif strcmp(comp1,'2')==1
    spec_1 = spec_mat(:,3);
    isgood1 = goodH;
elseif strcmp(comp1,'P')==1
    spec_1 = spec_mat(:,4);
    isgood1 = goodP;
end

if strcmp(comp2,'Z')==1
    spec_2 = spec_mat(:,1);
    isgood2 = goodZ;
elseif strcmp(comp2,'1')==1
    spec_2 = spec_mat(:,2);
    isgood2 = goodH;
elseif strcmp(comp2,'2')==1
    spec_2 = spec_mat(:,3);
    isgood2 = goodH;
elseif strcmp(comp2,'P')==1
    spec_2 = spec_mat(:,4);
    isgood2 = goodP;
end

if strcmp(comp3,'Z')==1
    spec_3 = spec_mat(:,1);
    isgood3 = goodZ;
elseif strcmp(comp3,'1')==1
    spec_3 = spec_mat(:,2);
    isgood3 = goodH;
elseif strcmp(comp3,'2')==1
    spec_3 = spec_mat(:,3);
    isgood3 = goodH;
elseif strcmp(comp3,'P')==1
    spec_3 = spec_mat(:,4);
    isgood3 = goodP;
end

if strcmp(comp4,'Z')==1
    spec_4 = spec_mat(:,1);
    isgood4 = goodZ;
elseif strcmp(comp4,'1')==1
    spec_4 = spec_mat(:,2);
    isgood4 = goodH;
elseif strcmp(comp4,'2')==1
    spec_4 = spec_mat(:,3);
    isgood4 = goodH;
elseif strcmp(comp4,'P')==1
    spec_4 = spec_mat(:,4);
    isgood4 = goodP;
end



corrspec1_2 = spec_1-((TFs(TFidx12).transfunc_tap)'.*spec_2);
% corrspec3_2 = spec_3-((TFs(TFidx32).transfunc_tap)'.*spec_1);
corrspec3_2 = spec_3-((TFs(TFidx32).transfunc_tap)'.*spec_2);
corrspec1_32 = corrspec1_2-((TFs(TFidx13_2).transfunc_tap)'.*corrspec3_2);
corrspec4_2 = spec_4-((TFs(TFidx42).transfunc_tap)'.*spec_2);
corrspec4_32 = corrspec4_2-((TFs(TFidx43_2).transfunc_tap)'.*corrspec3_2);
corrspec1_432 = corrspec1_32-((TFs(TFidx14_32).transfunc_tap)'.*corrspec4_32);

amp1_2 = real(ifft(2*corrspec1_2,NFFT)./dt);
amp1_32 = real(ifft(2*corrspec1_32,NFFT)./dt);
amp1_432 = real(ifft(2*corrspec1_432,NFFT)./dt);

corrtime1_2 = amp1_2(npad0+1:npad0+npts); 
corrtime1_32 = amp1_32(npad0+1:npad0+npts);
corrtime1_432 = amp1_432(npad0+1:npad0+npts);

taptime1=0.01;
[spec_Z,npad01,npts1,f1,NFFT1]=spectra(Zraw,dt,taxisZ,taptime1);
[spec_H1,npad01,npts1,f1,NFFT1]=spectra(H1raw,dt,taxisZ,taptime1);
[spec_H2,npad01,npts1,f1,NFFT1]=spectra(H2raw,dt,taxisZ,taptime1);
[spec_P,npad01,npts1,f1,NFFT1]=spectra(Praw,dt,taxisZ,taptime1);

Transfunc_tap1_2=interp1(f,(TFs(TFidx12).transfunc_tap),f1);
Transfunc_tap3_2=interp1(f,(TFs(TFidx32).transfunc_tap),f1);
Transfunc_tap13_2=interp1(f,(TFs(TFidx13_2).transfunc_tap),f1);
Transfunc_tap4_2=interp1(f,(TFs(TFidx42).transfunc_tap),f1);
Transfunc_tap43_2=interp1(f,(TFs(TFidx43_2).transfunc_tap),f1);
Transfunc_tap1_432=interp1(f,(TFs(TFidx14_32).transfunc_tap),f1);

 
corrspec1_2_day=spec_Z-(Transfunc_tap1_2'.*spec_H1);
corrspec3_2_day =spec_H2-(Transfunc_tap3_2'.*spec_H1);
corrspec1_32_day = corrspec1_2_day-(Transfunc_tap13_2'.*corrspec3_2_day);
corrspec4_2_day = spec_P-(Transfunc_tap4_2'.*spec_H1);
corrspec4_32_day = corrspec4_2_day-(Transfunc_tap43_2'.*corrspec3_2_day);
corrspec1_432_day = corrspec1_32_day-(Transfunc_tap1_432'.*corrspec4_32_day);
 
amp1_2_day = real(ifft(2*corrspec1_2_day,NFFT1)./dt);
amp1_32_day = real(ifft(2*corrspec1_32_day,NFFT1)./dt);
amp1_432_day = real(ifft(2*corrspec1_432_day,NFFT1)./dt);

corrtime1_2_day = amp1_2_day(npad01+1:npad01+npts1);
corrtime1_32_day = amp1_32_day(npad01+1:npad01+npts1);
corrtime1_432_day = amp1_432_day(npad01+1:npad01+npts1);
% interp(TFs(TFidx12).transfunc_tap)'
% T1=10;T2=150;fn = 1/2/dt;
% [b,a]=butter(2,[1/fn/T2,1/fn/T1]);
% figure;plot(filtfilt(b,a,Zraw),'b');
% hold on; plot(filtfilt(b,a,corrtime1_432_day),'r');
% 
% figure;plot(filtfilt(b,a,corrtime1_432_day(1:length(corrtime1_432))),'r');
% hold on;plot(filtfilt(b,a,corrtime1_432),'k');

noisespec1_2=((TFs(TFidx12).transfunc_tap)'.*spec_2);
noisespec3_2=((TFs(TFidx32).transfunc_tap)'.*spec_2);
noisespec4_2=((TFs(TFidx42).transfunc_tap)'.*spec_2);
% noisespec1_32=((TFs(TFidx13_2).transfunc_tap)'.*corrspec3_2);
% %noisespec1_32=((TFs(TFidx13_2).transfunc_tap)'.*spec_sig3_2);
% noisespec4_32=((TFs(TFidx43_2).transfunc_tap)'.*corrspec3_2);
% noisespec1_432=((TFs(TFidx14_32).transfunc_tap)'.*corrspec4_32);

ampnoise1_2=real(ifft(2*noisespec1_2,NFFT)./dt);
ampnoise3_2=real(ifft(2*noisespec3_2,NFFT)./dt);
ampnoise4_2=real(ifft(2*noisespec4_2,NFFT)./dt);
% ampnoise1_32=real(ifft(2*noisespec1_32,NFFT)./dt);
% ampnoise4_32=real(ifft(2*noisespec4_32,NFFT)./dt);
% ampnoise1_432=real(ifft(2*noisespec1_432,NFFT)./dt);

noise1_2=ampnoise1_2(npad0+1:npad0+npts); 
noise3_2=ampnoise3_2(npad0+1:npad0+npts);
noise4_2=ampnoise4_2(npad0+1:npad0+npts);
% noise1_32=ampnoise1_32(npad0+1:npad0+npts);
% noise4_32=ampnoise4_32(npad0+1:npad0+npts);
% noise1_432=ampnoise1_432(npad0+1:npad0+npts);

corrsingal1_2=Zsig-noise1_2;
corrsingal3_2=H2sig-noise3_2;
corrsignal4_2=Psig-noise4_2;
[spec_sig3_2,npad0,npts,f,NFFT]=spectra(corrsingal3_2,dt,taxisZ,taptime);
%[spec_sig4_2,npad0,npts,f,NFFT]=spectra(corrsignal4_2,dt,taxisZ,taptime);

noisespec1_32=((TFs(TFidx13_2).transfunc_tap)'.*spec_sig3_2);
noisespec4_32=((TFs(TFidx43_2).transfunc_tap)'.*spec_sig3_2);
ampnoise1_32=real(ifft(2*noisespec1_32,NFFT)./dt);
ampnoise4_32=real(ifft(2*noisespec4_32,NFFT)./dt);
noise1_32=ampnoise1_32(npad0+1:npad0+npts);
noise4_32=ampnoise4_32(npad0+1:npad0+npts);

%corrsingal1_32=corrsingal1_2-noise1_32;
corrsignal4_32=corrsignal4_2-noise4_32;

[spec_sig4_32,npad0,npts,f,NFFT]=spectra(corrsignal4_32,dt,taxisZ,taptime);

noisespec1_432=((TFs(TFidx14_32).transfunc_tap)'.*spec_sig4_32);
ampnoise1_432=real(ifft(2*noisespec1_432,NFFT)./dt);
noise1_432=ampnoise1_432(npad0+1:npad0+npts);

%corrsignal1_432=corrsingal1_32-noise1_432;


corrgood1_2 = isgood1*isgood2;
corrgood1_32 = isgood1*isgood2*isgood3;
corrgood1_432 = isgood1*isgood2*isgood3*isgood4;


return