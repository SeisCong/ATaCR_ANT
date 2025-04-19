function [spec,npad0,npts,f,NFFT]=spectra(signal,dt,taxis,taptime)

% Spectra parameters for resampled data
npts = length(signal);
Ppower = nextpow2(npts);
NFFT = 2^Ppower;
samprate = 1/dt; %common named Fs
f = samprate/2*linspace(0,1,NFFT/2+1);
pts_begin = 1;
pts_end = length(signal);
npad0 = floor((NFFT-npts)/2);

amp  = signal(pts_begin:pts_end);
amp  = amp.*flat_hanning(taxis(pts_begin:pts_end),taptime*dt*npts);
amp  = detrend(amp,0);
amp  = padarray(amp,[npad0 0],'both');
spectrum = fft(amp,NFFT).*dt;
spec = spectrum(1:NFFT/2+1);

return