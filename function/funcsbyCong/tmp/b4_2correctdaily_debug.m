%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correctdaily.m
% script for correcting daily record for tilt and compliance noise
% Modified from b4.correctevent.m from Helen
%% Cong Li
%% conli@geo.umass.edu
%% 07/09/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/Volumes/CongSeisDisk/MultiModeSurf_ANTomogrpahy/ATaCR-master/function/funcsbyCong')
addpath('/Volumes/CongSeisDisk/MultiModeSurf_ANTomogrpahy/ATaCR-master/function/funcsbyCong/tmp')
clear all; close all

isfigure_sta = 1;
isfigure_spectra = 1;
issavefigure = 1;
isoverwrite =1;

T1 = 10; T2= 150; %filter period range for plotting seismic data
Workfolder='/Volumes/CongSeisDisk/MultiModeSurf_ANTomogrpahy/Data/';

% input lists of stations with bad components, file doesn't have to exists
% leave as default if no bad stations
badstalist = 'NOISETC_CI/Bad_Z.txt';
badhorzlist = 'NOISETC_CI/Bad_H.txt';
badpreslist = 'NOISETC_CI/Bad_P.txt';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Loading parameters
setup_parameter;
% skipp station/channel list
if exist(badstalist)
	badstas = textread(badstalist,'%s');
else
    badstas = [];
end
if exist(badhorzlist)
	badhorz = textread(badhorzlist,'%s');
else
    badhorz = [];
end
if exist(badpreslist)
	badpres = textread(badpreslist,'%s');
else 
    badpres = [];
end
% Transfer list
for itf = 1:length(TF_list)
    TF_size(itf) = length(TF_list{itf});
end
[vec,TF_list_sort] = sort(TF_size,'descend');

inpath_dir = sprintf('%s/TRANSFUN/',OUTdir);
TF_stadir = dir(fullfile(inpath_dir));
ii=1;
% Load transfer function
for itf = 1:length(TF_stadir)
    station = TF_stadir(itf).name;
    inpath = sprintf('%s/TRANSFUN/%s/',OUTdir,station);
    if ~isdir(inpath)
        continue
    end
    filenames = dir(fullfile(inpath,['*.mat']));
    if isempty(filenames)
        continue
    end
    TF_stalist{ii} = station;
    ii=ii+1;
end
% Traverse each station folder
for istat=28%1:length(stations)
inpath_day = [Workfolder 'NOISETC_CI/DATA/datacache_day_preproc/YO/' stations{istat} '/']; % path to daily record
Day_filename=[];
Day_filename = dir(fullfile(inpath_day,['*.mat']));
for ie = 6:length(Day_filename) 
    close all
    clear corrseis_matrix
    % Load daily data
    if length(Day_filename(ie).name)~=23
        continue
    end
    eventid = Day_filename(ie).name(1:12);
    disp(['Day: ',eventid]);
    clear spec_mat
    TF_check = zeros(size(TF_size));
    idx1 = findstr(Day_filename(ie).name,'_');
    idx2 = findstr(Day_filename(ie).name,'.mat');     
    network = Day_filename(ie).name(idx1(1)+1:idx1(2)-1);
    station = Day_filename(ie).name(idx1(2)+1:idx2-1);   
    netsta = [network,station]; 
    TF_staidx = find(strcmp(netsta,TF_stalist)==1);
    if isempty(TF_staidx)
         continue
    end   
    sta = load(fullfile(inpath_day,Day_filename(ie).name)); %loading the correct file     
    % Specify transfer function file path and output paths
   inpath_trans = sprintf('%s/TRANSFUN/%s/',OUTdir,netsta);
    if tf_op ==1
       outpath = sprintf('%s/CORRSEIS/%s/',OUTdir,netsta);
       figoutpath=sprintf('%s/CORREVENTS',FIGdir);
    elseif tf_op ==2
       outpath = sprintf('%s/CORRSEISAVTF/%s/',OUTdir,netsta);
       figoutpath=sprintf('%s/CORREVENTS',FIGdir);
    end
		
   if ~isoverwrite && exist(sprintf('%s/%s_%s_corrseis.mat',outpath,netsta,eventid))==2
        disp([eventid,' exists. Skipping...'])
        continue
   end
   % Initialize the processed daily data           
   [Zraw,H1raw,H2raw,Praw,taxisZ,taxis1,taxis2,taxisP,dt] = varsetup_correctdailyrecord(...
       sta,chz_vec,ch1_vec,ch2_vec,chp_vec,T);
   if isnan(Zraw)
       disp('No Z component, skipping');
       continue;
   end     
   % Finding and loading the tranfer files; 
   % Note: tf_op=2 average transfer file for all data recorded at a station
   %       tf_op=1 daily transfer file
    trans_filename=correctdaily_findTFs(inpath_trans,tf_op,netsta,eventid);
     if isnan(trans_filename)
       disp('No transfer function file, skipping')
         continue
    end    
        if exist(trans_filename,'file')
            goodP = 1;
            goodH=1;
            goodZ=1;
            if ~isempty(badstas)
            if ~isempty(find(~cellfun('isempty', strfind(badstas,netsta)))) == 1
                disp('Vertical Bad.')
                goodZ=0; %if pressure is marked as bad, turn good flag off
            end
            end
            if ~isempty(badpres)
            if ~isempty(find(~cellfun('isempty', strfind(badpres,netsta)))) == 1
                disp('Pressure Bad.')
                goodP=0; %if pressure is marked as bad, turn good flag off
            end
            end
            if ~isempty(badhorz)
            if ~isempty(find(~cellfun('isempty', strfind(badhorz,netsta)))) == 1
                disp('Horizontal Bad.')
                goodH=0; %if horizontal is marked as bad, turn good flag off
            end
            end
        end
      disp('Calculating Corrected Seismograms...');
      if tf_op==1
      disp(sprintf('Using TF: %s', trans_filename(end-30:end)));
      elseif tf_op==2
      disp(sprintf('Using TF: %s', trans_filename(end-25:end)));
      end
%% Start to remove tilt and/or compliance noise for each daily file 
       taperwin=0;%round(T/(dt*4));
       sm_factor=1500; win_s=500;win_ss=200;
 for iwindow=1:round((length(Zraw)-(T/dt))/((T/dt)-taperwin))+1
     %length(Zraw)/(T/dt)
% Calculate event spectra properties for each window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    debug
            fn = 1/2/dt;
            [b,a]=butter(2,[1/fn/T2,1/fn/T1]);
        if iwindow==2
            ab=(T/dt)*(iwindow-2)+1:(T/dt)*(iwindow-1);
            signal=[];
            signal=corrseis(1).timeseries;
            figure; plot(ab,filtfilt(b,a,signal),'r');
          residual=Zraw(1:length(corrseis(1).noise))-corrseis(1).noise;
          %-corrseis(3).noise-corrseis(2).noise
           hold on;plot(ab,filtfilt(b,a,residual),'b');
        elseif iwindow>2 
           ab=1:(T/dt)*(iwindow-1)-(iwindow-2)*taperwin;
           %(T/dt)*(iwindow-2)+1-(iwindow-2)*taperwin:...
            %   (T/dt)*(iwindow-1)-(iwindow-2)*taperwin;
          signal=[];
          signal=corrseis(3).timeseries;
           hold on;plot(ab,filtfilt(b,a,signal));
          %residual=Zraw(1:length(corrseis(1).noise))-corrseis(3).noise-corrseis(2).noise-corrseis(1).noise;
 
          % hold on;plot(ab,filtfilt(b,a,residual));
           hold on;plot(ab,filtfilt(b,a,corrtime1_432_day(1:length(corrseis(3).noise))));
          % figure;plot(ab,filtfilt(b,a,Zraw(1:length(corrseis(1).noise))));
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if iwindow==1
        Zsig=Zraw((T/dt)*(iwindow-1)+1:(T/dt)*iwindow);
        H1sig=H1raw((T/dt)*(iwindow-1)+1:(T/dt)*iwindow);
        H2sig=H2raw((T/dt)*(iwindow-1)+1:(T/dt)*iwindow);
        Psig=Praw((T/dt)*(iwindow-1)+1:(T/dt)*iwindow);
        elseif iwindow==round((length(Zraw)-(T/dt))/((T/dt)-taperwin))+1;%length(length(Zraw)/(T/dt))
        Zsig=Zraw(end-(T/dt)+1:end);
        H1sig=H1raw(end-(T/dt)+1:end);
        H2sig=H2raw(end-(T/dt)+1:end);
        Psig=Praw(end-(T/dt)+1:end);
        else
        Zsig=Zraw((T/dt)*(iwindow-1)+1-(iwindow-1)*taperwin:...
            (T/dt)*iwindow-(iwindow-1)*taperwin); 
        H1sig=H1raw((T/dt)*(iwindow-1)+1-(iwindow-1)*taperwin:...
            (T/dt)*iwindow-(iwindow-1)*taperwin); 
        H2sig=H2raw((T/dt)*(iwindow-1)+1-(iwindow-1)*taperwin:...
            (T/dt)*iwindow-(iwindow-1)*taperwin); 
        Psig=Praw((T/dt)*(iwindow-1)+1-(iwindow-1)*taperwin:...
            (T/dt)*iwindow-(iwindow-1)*taperwin); 
        end
        % Analysis spectrum feature of a signal (Z,H1,H2,P)
        [spec_Z,npad0,npts,f,NFFT]=spectra(Zraw,dt,taxisZ,taptime);
        [spec_mat,npad0,npts,f,NFFT]=spectra_correctevent(Zsig,...
            H1sig,H2sig,Psig,dt,taxisZ,taxis1,taxis2,taxisP,taptime);      
        if ~exist(figoutpath)
            mkdir(figoutpath);
        end
        if ~exist(outpath)
            mkdir(outpath);
        end 
        % Load transfer functions
        load(trans_filename);
        freqcomp=transprop.params.freqcomp;
        if filop==1
            lp = taper_lim(2);
            hp = taper_lim(1);
        elseif filop ==2
            lp = freqcomp+0.005;
            hp=0;
        end
        NFFT = transprop.params.NFFT;
        f=transprop.params.f;
        for itf = 1:length(TFs); %taper each of the individual transfer functions
            TFs(itf).transfunc_tap = tf_taper((TFs(itf).transfunc)',f',hp,lp,tap_width);
        end
% Start to remove compliance and/or tilt noise using each combination of transfer function 
        TF_check = zeros(size(TF_size)); %% L.C.
        ii = 1;
        for itf = 1:length(TF_list_sort)
            TF_cal = TF_list(TF_list_sort(itf));
            if TF_check(TF_list_sort(itf)) == 0 % check if TF calculated already
                % check if rotational or not
                if ~isempty(strfind(cell2mat(TF_cal),'H'))
                if  tf_op == 2 % average doesn't make sense for this
                    disp('Average transfer function is not appropriate for removal noise with H')
                    continue
                end
                if length(TF_cal{1}) == 2 % 1 component rotational TF, i.e. ZH
                    [corrspec,corrtime,corrgood,label_list] = tfcomp1rotate_correctevent(TF_cal,TFs,transprop,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);
                    TF_check(TF_list_sort(itf)) = 1;
                    corrseis(ii).label = label_list{1};
                    if ~isfield(corrseis(ii),'spectrum')  
                    corrseis(ii).spectrum = corrspec;
                    else
                    corrseis(ii).spectrum =[corrseis(ii).spectrum, corrspec];
                    end
                    if ~isfield(corrseis(ii),'timeseries')
                    corrseis(ii).timeseries = corrtime;
                    else
                    corrseis(ii).timeseries = ...
                        merge_dailycorreis(corrseis(ii).timeseries, corrtime,taperwin);    
                    end
                    corrseis(ii).isgood = corrgood;
                    ii = ii+1;
                elseif length(TF_cal{1}) == 4 % ZP-H
                    [corrspec1_2,corrspec1_32,corrtime1_2,corrtime1_32,corrgood1_2,corrgood1_32,...
                        label_list] = tfcomp2rotate_correctevent(TF_cal,TFs,transprop,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);
                    TF_check(TF_list_sort(itf)) = 1;
                    corrseis(ii).label = label_list{1};
                if  ~isfield(corrseis(ii),'spectrum')  
                    corrseis(ii).spectrum = corrspec1_2;
                else
                    corrseis(ii).spectrum = [corrseis(ii).spectrum, corrspec1_2];
                end
                if ~isfield(corrseis(ii),'timeseries') | isempty(corrseis(ii).timeseries)
                    corrseis(ii).timeseries = corrtime1_2;
                else
                    corrseis(ii).timeseries = merge_dailycorreis(corrseis(ii).timeseries,corrtime1_2,taperwin);
                end
                    corrseis(ii).isgood = corrgood1_2;
                    if ~isempty(find(strcmp(label_list{1},TF_list)==1))
                        TFidx = find(strcmp(label_list{1},TF_list)==1);
                        TF_check(TFidx)=1;
                    end
                    corrseis(ii+1).label = label_list{2};
                    if  ~isfield(corrseis(ii+1),'spectrum') 
                        corrseis(ii+1).spectrum = corrspec1_32;
                    else
                        corrseis(ii+1).spectrum = [corrseis(ii+1).spectrum,corrspec1_32];
                    end
                    if ~isfield(corrseis(ii+1),'timeseries') | isempty(corrseis(ii+1).timeseries)
                    corrseis(ii+1).timeseries = corrtime1_32;
                    else
                    corrseis(ii+1).timeseries =...
                        merge_dailycorreis(corrseis(ii+1).timeseries, corrtime1_32,taperwin);
                    end
                    corrseis(ii+1).isgood = corrgood1_32;
                    if ~isempty(find(strcmp(label_list{2},TF_list)==1))
                        TFidx = find(strcmp(label_list{2},TF_list)==1);
                        TF_check(TFidx)=1;
                    end
                    ii = ii+2;
                end
                else % if component wise (No rotation of H)
                if length(TF_cal{1}) == 2 % 1 component correction
                    [corrspec,corrtime,corrgood,label_list] = tfcomp1_correctevent(TF_cal,TFs,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);
                    TF_check(TF_list_sort(itf)) = 1;
                    corrseis(ii).label = label_list{1};
                if ~isfield(corrseis(ii),'spectrum')    
                    corrseis(ii).spectrum = corrspec;
                else
                    corrseis(ii).spectrum = [corrseis(ii).spectrum,corrspec];
                end
                if ~isfield(corrseis(ii),'timeseries') | isempty(corrseis(ii).timeseries)
                    corrseis(ii).timeseries = corrtime;
                else
                    corrseis(ii).timeseries = ...
                        merge_dailycorreis(corrseis(ii).timeseries, corrtime,taperwin);
                   % [corrseis(ii).timeseries; corrtime];
                end           
                    corrseis(ii).isgood = corrgood;
                    ii = ii+1;
                elseif length(TF_cal{1}) == 5 % 3 component correction
%                     [corrspec1_2,corrspec1_32,corrspec1_432,corrtime1_2,...
%                         corrtime1_32,corrtime1_432,corrgood1_2,corrgood1_32,corrgood1_432,label_list] =...
%                         tfcomp3_correctevent(TF_cal,TFs,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);
                    
                    [corrspec1_2,corrspec1_32,corrspec1_432,corrtime1_2,corrtime1_32,corrtime1_432,...
                        corrtime1_2_day,corrtime1_32_day,corrtime1_432_day,corrgood1_2,corrgood1_32,corrgood1_432,label_list,...
                    noise1_2,noise3_2,noise1_32,noise4_32,noise1_432] =...
                    tfcomp3_correctdaily(TF_cal,TFs,spec_mat,Zraw,H1raw,H2raw,Praw,Zsig,H1sig,H2sig,Psig,...
                    NFFT,dt,f,npad0,npts,goodP,goodH,goodZ,taxisZ,taptime);
         
                    %Z-H1
                    TF_check(TF_list_sort(itf)) = 1;% 1 means already done
                    corrseis(ii).label = label_list{1};
                    if ~isfield(corrseis(ii),'spectrum')
                    corrseis(ii).spectrum = corrspec1_2;  
                    else
                    corrseis(ii).spectrum = [corrseis(ii).spectrum,corrspec1_2];
                    end
                    
                    if ~isfield(corrseis(ii),'timeseries') | isempty(corrseis(ii).timeseries)
                    corrseis(ii).timeseries = corrtime1_2;
                    else
                    corrseis(ii).timeseries = ...
                        merge_dailycorreis(corrseis(ii).timeseries,corrtime1_2,taperwin);
                    end
                    
                    if ~isfield(corrseis(ii),'noise')
                    corrseis(ii).noise = noise1_2;
                    else
                    corrseis(ii).noise = merge_dailycorreis(corrseis(ii).noise,noise1_2,taperwin);
                    end
                    %% SMOTHING PART L.C.
                    if iwindow>=2
                    sm_win=[];int_win=[];tt_win=[];
                    tt_win1=[(iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)-win_ss];
                    tt_win2=[(iwindow-1)*(T/dt)+win_ss:(iwindow-1)*(T/dt)+win_s];
                    tt_win=[tt_win1,tt_win2];
                    tt_win1=(iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s;
                    int_win=interp1(tt_win,corrseis(ii).noise(tt_win),tt_win1,'cubic');
                    sm_win=smooth(corrseis(ii).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s),sm_factor,'loess');
                    figure;plot(corrseis(ii).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s));
                    corrseis(ii).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s)=int_win;%sm_win;
                    hold on;plot(corrseis(ii).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s));
                    close all
                    end
                    %%
                    corrseis(ii).isgood = corrgood1_2;
                    
                    % check and declear that removal has been done 
                    if ~isempty(find(strcmp(label_list{1},TF_list)==1))
                        TFidx = find(strcmp(label_list{1},TF_list)==1); 
                        TF_check(TFidx)=1;
                    end
                    %Z-H1-H2
                    corrseis(ii+1).label = label_list{2};
                    if ~isfield(corrseis(ii+1),'spectrum')
                    corrseis(ii+1).spectrum = corrspec1_32;
                    else
                    corrseis(ii+1).spectrum =[corrseis(ii+1).spectrum,corrspec1_32];
                    end
                    
                    if ~isfield(corrseis(ii+1),'timeseries') | isempty(corrseis(ii+1).timeseries)
                    corrseis(ii+1).timeseries = corrtime1_32;    
                    else
                    corrseis(ii+1).timeseries = ...
                        merge_dailycorreis(corrseis(ii+1).timeseries,corrtime1_32,taperwin);
                    %[corrseis(ii+1).timeseries;corrtime1_32];
                    end
                    
                     if ~isfield(corrseis(ii+1),'noise')
                    corrseis(ii+1).noise = noise1_32;
                    else
                    corrseis(ii+1).noise = merge_dailycorreis(corrseis(ii+1).noise,noise1_32,taperwin);
                    end
                    %% SMOTHING PART L.C.
                    if iwindow>=2
                    sm_win=[];int_win=[];tt_win=[];
                    tt_win1=[(iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)-win_ss];
                    tt_win2=[(iwindow-1)*(T/dt)+win_ss:(iwindow-1)*(T/dt)+win_s];
                    tt_win=[tt_win1,tt_win2];
                    tt_win1=(iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s;
                    int_win=interp1(tt_win,corrseis(ii+1).noise(tt_win),tt_win1,'cubic');
           %         sm_win=smooth(corrseis(ii+1).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s),sm_factor,'loess');
                    figure;plot(corrseis(ii+1).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s));
                    corrseis(ii+1).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s)=int_win;%sm_win;
                    hold on;plot(corrseis(ii+1).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s));
                    close all
                    end
                    
                    corrseis(ii+1).isgood = corrgood1_32;
                    if ~isempty(find(strcmp(label_list{2},TF_list)==1))
                        TFidx = find(strcmp(label_list{2},TF_list)==1);
                        TF_check(TFidx)=1;
                    end
                    %Z-H1-H2-P
                    corrseis(ii+2).label = label_list{3};
                    if  ~isfield(corrseis(ii+2),'spctrum')
                    corrseis(ii+2).spectrum = corrspec1_432;
                    else
                    corrseis(ii+2).spectrum =[corrseis(ii+2).spectrum,corrspec1_432];
                    end
                    if ~isfield(corrseis(ii+2),'timeseries') | isempty(corrseis(ii+2).timeseries)
                    corrseis(ii+2).timeseries=corrtime1_432;
                    else
                    corrseis(ii+2).timeseries = ...
                        merge_dailycorreis(corrseis(ii+2).timeseries,corrtime1_432,taperwin);
                    end
                    
                      if ~isfield(corrseis(ii+1),'noise')
                    corrseis(ii+2).noise = noise1_432;
                    else
                    corrseis(ii+2).noise = merge_dailycorreis(corrseis(ii+2).noise,noise1_432,taperwin);
                    end
                    %% SMOTHING PART L.C.
                    if iwindow>=2
                    sm_win=[];int_win=[];tt_win=[];
                    tt_win1=[(iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)-win_ss];
                    tt_win2=[(iwindow-1)*(T/dt)+win_ss:(iwindow-1)*(T/dt)+win_s];
                    tt_win=[tt_win1,tt_win2];
                    tt_win1=(iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s;
                    int_win=interp1(tt_win,corrseis(ii+2).noise(tt_win),tt_win1,'cubic');
           %         sm_win=smooth(corrseis(ii+1).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s),sm_factor,'loess');
                    figure;plot(corrseis(ii+2).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s));
                    corrseis(ii+2).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s)=int_win;%sm_win;
                    hold on;plot(corrseis(ii+2).noise((iwindow-1)*(T/dt)-win_s:(iwindow-1)*(T/dt)+win_s));
                    close all
                    end
                    
                    corrseis(ii+2).isgood = corrgood1_432;
                    if ~isempty(find(strcmp(label_list{3},TF_list)==1))
                        TFidx = find(strcmp(label_list{3},TF_list)==1);
                        TF_check(TFidx)=1;
                    end
                    ii=ii+3;
                end
                end
            end
        end
 end
    
 
        corrected.params.f = f;
        corrected.params.tf_op = tf_op;
        corrected.params.eventid = eventid;
        corrected.params.filop = filop;
        corrected.params.taxis = taxisZ;
        corrected.params.station = station;
        corrected.params.network = network;
        corrected.params.NFFT = NFFT;
        corrected.params.dt = dt;
        corrected.params.TFfilename = trans_filename;
%         corrected.params.stalat = 
%         corrected.params.stalon = 
%         corrected.params.freqcomp = 
        
        filename = sprintf('%s/%s_%s_corrseis.mat',outpath,netsta,eventid);
        save(filename,'corrseis','corrected');
        
        %station plots
        if isfigure_sta == 1;
            spitplots_correctevent(dt,T1,T2,Zraw,H1raw,H2raw,Praw,taxisZ,eventid,netsta,corrseis)
            if issavefigure==1
                figure(101)
                filename=sprintf('%s/%s_%s_originalseis',figoutpath,eventid,netsta);
                print(gcf,'-dpng',filename)
                if tf_op ==1
                    figure(102)
                    filename=sprintf('%s/%s_%s_corrseis',figoutpath,eventid,netsta);
                    print(gcf,'-dpng',filename)
                elseif tf_op==2
                    figure(102)
                    filename=sprintf('%s/%s_%s_corrseis_av',figoutpath,eventid,netsta);
                    print(gcf,'-dpng',filename)
                end
            end
        end
		%plot spectra against Peterson's low and high noise models
		if isfigure_spectra
            plot_corrspectra(T1,T2,Zraw,corrseis,f,NFFT,dt,taxisZ,eventid,netsta,lp,hp)
            if issavefigure==1
                if tf_op ==1
                    figure(1)
                    filename=sprintf('%s/%s_%s_spectra',figoutpath,eventid,netsta);
                    print(gcf,'-dpng',filename)
                elseif tf_op==2
                    figure(1)
                    filename=sprintf('%s/%s_%s_spectra_av',figoutpath,eventid,netsta);
                    print(gcf,'-dpng',filename)
                end
            end
        end
   % end
    
end
end
