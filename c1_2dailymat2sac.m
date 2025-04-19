% Corrected daily records into sac format
%
% Convert *.mat files containing corrected seismograms to SAC files matching
% the original input SAC files.
%
% Cong Li
% conli@geo.umass.edu
% updated 07/21

clear all;clc;
addpath('/home/licong/Project/Software/Mat_Lib/Common')
setup_parameter;

inpath_uncorr = '/media/licong/EastTibet/Data_ENAM/Data/YO_FTN_old'; % Originial sac files
str_corr = 'ZP-21'; % String for correction to output
channel = 'BHZ';
fid=fopen([OUTdir '/CORRSEIS_SAC/MissConvert.txt'],'w');
%% Load data
if tf_op == 1
    corrseis_path = sprintf('%s/CORRSEIS/',OUTdir);
elseif tf_op ==2
    corrseis_path = sprintf('%s/CORRSEISAVTF/',OUTdir);
end

stadirs = dir(fullfile(corrseis_path));
for ista = 1:length(stadirs)
    station = stadirs(ista).name;
    inpath_corr = sprintf('%s%s/',corrseis_path,station);
    if length(station)~=5%~isdir(inpath_corr)
        continue
    end
    filenames_corr = dir(fullfile(inpath_corr,['*.mat']));
    disp(station);
    % Loop over event files
    for iev = 1:length(filenames_corr);
        if ~exist(fullfile(inpath_corr,filenames_corr(iev).name))
            continue
        end
        load(fullfile(inpath_corr,filenames_corr(iev).name))
        corr_idx = find(strcmp({corrseis.label},str_corr));
        corrdata = corrseis(corr_idx).timeseries;
        
        
        yearcorr=corrected.params.eventid(1:4);
        monthcorr=corrected.params.eventid(5:6);
        datecorr=corrected.params.eventid(7:8);
        hourSec=corrected.params.eventid(9:end);
        
        daycorr=datenum([yearcorr '-' monthcorr '-' datecorr],'yyyy-mm-dd');
        daycorr_s=datenum([yearcorr '-1-1'],'yyyy-mm-dd');
        % Load data headers
        sacfilename=fullfile(sprintf('%s/%s.%s/ftn.%s.%s.%s.%s.HHZ.SAC',...
            inpath_uncorr,corrected.params.network, corrected.params.station, ...
            yearcorr, num2str(daycorr-daycorr_s+1,'%03.f'), corrected.params.network, corrected.params.station));
        if exist(sacfilename,'file')
        sacin = rdsac(sacfilename);
        if tf_op ==1
            opath = sprintf('%s/CORRSEIS_SAC/%s.%s/',OUTdir,corrected.params.network,...
                corrected.params.station);
         elseif tf_op ==2
            opath = sprintf('%s/CORRSEISAVTF_SAC/%s.%s/',OUTdir,corrected.params.network,...
                corrected.params.station);
        end
        if ~exist(opath)
            mkdir(opath);
        end
        H = sacin.HEADER;
        H.DELTA = corrected.params.dt;
        data1 = corrseis(corr_idx).timeseries;
        if length(corrected.params.taxis)<86400/corrected.params.dt+1
            H.NPTS=86400/corrected.params.dt+1; % fill zeros
            data = [data1;zeros(H.NPTS-length(corrected.params.taxis),1)];
        elseif length(corrected.params.taxis)>86400/corrected.params.dt+1
            H.NPTS=86400/corrected.params.dt+1; % fill zeros
            data = data1(1:H.NPTS);
        else
            H.NPTS = length(corrected.params.taxis);
        end
        clear data1
        evid = corrected.params.eventid;
        startDate = datetime(H.NZYEAR,1,H.NZJDAY); % use this to convert jday to month - day
        daysac=datenum(datestr(startDate),'dd-mm-yyyy');
        if daysac~=daycorr
              disp('Warining:NZJDAY is not equal! Correct the absolute time now')
                 H.NZJDAY=daycorr-daycorr_s+1;
                 H.NZHOUR=0;
                 H.NZMIN=0;
                 H.NZSEC=0;
                 H.NZMSEC=0;
                 startDate = datetime(H.NZYEAR,1,H.NZJDAY);
        end
        if H.NZMSEC~=0 
            H.NZMSEC=0;
        end
        startYear = num2str(year(startDate), '%04d');
        startMonth = num2str(month(startDate), '%02d');
        startDay = num2str(day(startDate), '%02d');
        
      fullevid = [startYear,startMonth,startDay,num2str(H.NZHOUR,'%02d'),...
            num2str(H.NZMIN,'%02d'),num2str(H.NZSEC,'%02d'),num2str(H.NZMSEC,'%03d')];
        startTime = datenum(fullevid,'yyyymmddhhMMSSFFF');
        sac_path = fullfile(sprintf('%s/%s.%s.%s.%s.%s.SAC',opath,yearcorr,...
            num2str(daycorr-daycorr_s+1,'%03.f'), sacin.HEADER.KNETWK, sacin.HEADER.KSTNM,sacin.HEADER.KCMPNM));
        mksac(sac_path,data,startTime,H); 
        disp(startDate)
        %disp([num2str(H.NZYEAR) '-' num2str(H.NZJDAY)])
        else 
            fprintf(fid, 'did not find file %s\n',sacfilename);
        end
    end
 end
fclose(fid)