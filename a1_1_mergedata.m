clc;clear all;close all;
addpath('/media/licong/EastTibet/Data_ENAM/ATaCR-master/function/funcsbyCong')
Workfolder='/media/licong/EastTibet/Data_ENAM/';
INPUTdir = [Workfolder 'NOISETC_CI/DATA/datacache_day/'];
network='YO';
network1='YO_ZP';
network2='YO_HH';
station = textread([Workfolder 'NOISETC_CI/Stationlist.txt'],'%s'); 
%station='X05';
fid=fopen('StationForMerge.txt','w');
for i=1:length(station)
filesuff = sprintf('*_%s_%s.mat',network,station{i});
data_filenames = dir(fullfile(INPUTdir,network1,station{i},'/',filesuff));
for ie =1:length(data_filenames)
disp(data_filenames(ie).name);
if ~exist([INPUTdir network1 '/' station{i} '/',data_filenames(ie).name],'file')
    disp('No vertical and pressure components! Skipping!');  
    continue;
else
    clear traces_day_ZP
    traces_day_ZP=load(fullfile(INPUTdir,network1,station{i},...
       '/',data_filenames(ie).name));
 channels={};
    for j=1:length(traces_day_ZP.traces_day)
    channels{j}=traces_day_ZP.traces_day(j).channel;
    end
    BDH_idx=find(strcmp(channels,'BDH')==1);
    HHZ_idx=find(strcmp(channels,'HHZ')==1);
    traces_day_ZP1=[];
    %% BDH merge
    if length(BDH_idx)==1
        traces_day_ZP1.traces_day(1)=traces_day_ZP.traces_day(BDH_idx);
    elseif length(BDH_idx)==0
        disp('No pressure component! Skipping');
        continue;
    else
        disp([num2str(length(BDH_idx)) ...
            ' segments in BDH component, will be merged now!']);    
    fprintf(fid,'station for merge: %s %s\n',station{i},data_filenames(ie).name);
       sep_data=[];
       sep_data=traces_day_ZP.traces_day(BDH_idx);
       traces_day_ZP1.traces_day=MergeSAC(sep_data);
       if isempty(traces_day_ZP1(1).traces_day)
         continue;
       end
    end
    %% HHZ merge   
    if length(HHZ_idx)==1
    elseif length(HHZ_idx)==0
        disp('No Z component! Skipping');
        continue;
    else
        disp([num2str(length(HHZ_idx)) ...
            ' segments in HHZ component, will be merged now!']); 
       sep_data=[];TEMP=[];
       sep_data=traces_day_ZP.traces_day(HHZ_idx);
       TEMP=MergeSAC(sep_data);
       if isempty(TEMP)
         continue;
       end
    traces_day_ZP1.traces_day=[traces_day_ZP1.traces_day,TEMP];
        clear traces_day_ZP
        traces_day_ZP=traces_day_ZP1;
        clear traces_day_ZP1
    end
end

if ~exist([INPUTdir network2 '/' station{i} '/',data_filenames(ie).name],'file')
    disp('No horizontal components!');
    continue;
else
clear traces_day_HH;
traces_day_HH=load(fullfile(INPUTdir,network2,station{i},'/',data_filenames(ie).name));    
channels={};
for j=1:length(traces_day_HH.traces_day)
    channels{j}=traces_day_HH.traces_day(j).channel;
end
 HH1_idx=find(strcmp(channels,'HH1')==1);
 HH2_idx=find(strcmp(channels,'HH2')==1);
 traces_day_HH1=[];
 %% HH1 merge
 if length(HH1_idx)==1
    elseif length(HH1_idx)==0
        disp('No pressure component! Skipping');
        continue;
 else
     
 disp([num2str(length(HH1_idx)) ...
       ' segments in HH1 component, will be merged now!']); 
       sep_data=[];
       sep_data=traces_day_HH.traces_day(HH1_idx);
       traces_day_HH1.traces_day=MergeSAC(sep_data);
       if isempty(traces_day_HH1(1).traces_day)
         continue;
       end
 end
 %% H2 merge
 if length(HH2_idx)==1
    elseif length(HH2_idx)==0
        disp('No pressure component! Skipping');
        continue;
 else
 disp([num2str(length(HH2_idx)) ...
       ' segments in HH2 component, will be merged now!']);    
       sep_data=[];TEMP=[];
       sep_data=traces_day_HH.traces_day(HH2_idx);
       TEMP=MergeSAC(sep_data);
       traces_day_HH1.traces_day=[traces_day_HH1.traces_day,TEMP];
       if isempty(TEMP)
         continue;
       end
        clear traces_day_HH
        traces_day_HH=traces_day_HH1;
        clear traces_day_HH1 
    end
end

clear traces_day;
traces_day=[traces_day_ZP.traces_day,traces_day_HH.traces_day];

if ~exist([INPUTdir network '/' station{i}] ,'dir')
    mkdir([INPUTdir network '/' station{i}])
end
save([INPUTdir network '/' station{i} '/' data_filenames(ie).name],'traces_day');


end
end
fclose(fid);