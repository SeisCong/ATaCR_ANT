% parameters for the NoiseTiltComp_pck

% path for matlab codes and functions
addpath ('function');
% location of the continuous matlab files for spectral properties
Workfolder='/media/licong/EastTibet/Data_ENAM/';
WORKINGdir = [Workfolder 'NOISETC_CI/DATA/datacache_day_preproc/'];
% output directory for spectra
OUTdir = [Workfolder 'NOISETC_CI/DATA/NOISETC'];
 % directory for figure output
FIGdir = [Workfolder 'NOISETC_CI/FIGURES/NOISETC'];


% information for station
network = 'YO';
stations = textread([Workfolder '/NOISETC_CI/Stationlist.txt'],'%s');

% Channel Names
chz_vec = {['HHZ'], ['BHZ'], ['LHZ']}; % list of acceptable names for Z component
ch1_vec = {['HH1'], ['BH1'], ['LHN'], ['LH1']}; % list of acceptable names for H1 component
ch2_vec = {['HH2'], ['BH2'], ['LHE'], ['LH2']}; % list of acceptable names for H2 component
chp_vec = {['HDH'], ['BDH']}; % list of acceptable names for P component

allchans = 1; % flag for if all channels are needed, or will process with missing channels (in trial period)

% Spectral Properties Windowing
T    = 7200;  % the legnth of each time window, in sec, should match the event data length for corrections
overlap = 0.3; %fraction of window overlap for spectra calculation

% Quality Control Parameters for Daily Windows
pb = [0.004 .2]; %pass-band, in Hz
tolerance = 1.5;%1.5;
a_val = 0.05;
minwin = 10;    % minimum numbers of time window for segment to be accepted

% Tilt orientation - only matters if using transfer functions with the 'H'
% option, but package needs variables specified to run; leave as default if
% not using
tiltfreq = [.005, .035]; % specifying frequency ranges for maximum coherence search
%tiltfreq = [1/20, 1/5]; % specifying frequency ranges for maximum coherence search

% Quality Control Parameters for Deployment Days
pb_dep = [0.004 .2]; %pass-band, in Hz
tolerance_dep = 3;%2;
a_val_dep = 0.05;

% Transfer Function Options
TF_list = {'ZP-21'};%{'ZP','Z1','Z2-1','ZP-21','ZH','ZP-H'};
% TF_list = {'ZP','Z1-P','Z2-1P','ZH','ZP-H'};
% convention: Z1 = transfer function between Z and H1
%             Z2-1 = transfer function between Z with H1 removed, and H2
%		      ZP-21 = transfer function between Z with H1 and H2 removed, and P
%             ZH = transfer function between Z and rotated max horizontal noise direction (H)

% Correction Options
taptime = 0.01;%0.075; % taper for seismogram correction step
                 % taperwin=taptime*dt*npts

tf_op = 1; %option for using either daily (1) or station average (2) TF in correction

filop = 1; %how to filter TF
% 1 - user specified constant high pass and low pass
% 2 - %lowpass - 0.005+freqcomp, adaptive to the infragravity wave cutoff, no high pass;
tap_width = 0.01; %width in percent of frequency vector of cosine taper
taper_lim = [0 0.6]; % upper and lower freuqncy in Hz limit for taper if user specified (option 1); zero means not applied
