function timeseriesMerge = merge_dailycorreis(timeseriesA,timeseriesB,taperwin)
% merge A and B, replace the overlaps part with the mean
%% Cong Li
%% conli@geo.umass.edu
%% 07/10/2021
timeseriesMerge=zeros(length(timeseriesA)+length(timeseriesB)-taperwin,1);
smoothwin=100;
if taperwin>length(timeseriesA)
    taperwin=length(timeseriesA);   
end
if taperwin==0
timeseriesMerge=[timeseriesA;timeseriesB];
else
tapercut=500; 
Meanpart=(timeseriesA(end-taperwin+1:end)+timeseriesB(1:taperwin))/2;
% AA=Meanpart;
%Meanpart=smooth(Meanpart,smoothwin);
% %avoid cut effect
% AB=Meanpart;
Meanpart(1:tapercut)=timeseriesA(end-taperwin+1:end-taperwin+tapercut);
Meanpart(end-tapercut+1:end)=timeseriesB(taperwin-tapercut+1:taperwin);

timeseriesMerge(1:length(timeseriesA)-taperwin)=...
    timeseriesA(1:length(timeseriesA)-taperwin);
timeseriesMerge(length(timeseriesA)-taperwin+1:length(timeseriesA))=Meanpart;
timeseriesMerge(length(timeseriesA)+1:end)=timeseriesB(taperwin+1:end);
end
%timeseriesMerge_filt=smooth(timeseriesMerge,smoothwin);
% figure;plot(timeseriesA);
% hold on;plot(7201:(7200+14400), timeseriesB);
% hold on;plot(timeseriesMerge);
return
end

