% make_starttimes
%
% Make start time files for daily data for ambient noise in 
% format YYYYmmddhhMM
%
% Cong Li; conli@geo.umass.edu
% updated 04/28/21
clc;clear all;
Workfolder='Project/MultiModeSurf_ANTomogrpahy/Data/';
if (exist([Workfolder 'NOISETC_CI'],'dir')~=1)
  mkdir([Workfolder 'NOISETC_CI']);
end
NetWork='YO'; %Network
starttime='201401010000'; % Start time for continuous data
endtime='201512310000'; % End time for continuous data
starttime_id=datenum(starttime,'yyyymmddHHMM');
endtime_id=datenum(endtime,'yyyymmddHHMM');
Duration=endtime_id-starttime_id+1;

for iday = 1:Duration
        day = datestr(starttime_id+iday-1,'yyyymmddHHMM'); % Assume files start at beginning of day
        daylist{iday} = day;
end

% Save starttimes to text file
dayfile=[Workfolder 'NOISETC_CI/' NetWork '_starttimes_CItest.txt'];
fid = fopen(dayfile,'w');
for iday = 1:length(daylist)
    fprintf(fid,'%s\n',daylist{iday});
end
fclose(fid);