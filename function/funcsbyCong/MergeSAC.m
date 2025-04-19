function temp_merge=MergeSAC(sep_data)
%% Merge seismic data into a daily record
       tick_skip=0;
       for h=1:length(sep_data)-1
          if ~isequal(sep_data(h).sacpz,sep_data(h+1).sacpz)
              disp('The PZ file varies, skipping!');
              return;
          end
       end
       for h=1:length(sep_data)-1
          if ~isequal(sep_data(h).sampleRate,sep_data(h+1).sampleRate)
              disp('sample rate varies, skipping!');
              return;
          end
       end
       temp_merge=sep_data(1); 
%        temp_merge.data=sep_data(1);
       temp_STime=datestr(sep_data(1).startTime,'yyyymmddHHMMss');
       STime=str2num(temp_STime(9:10))*3600+str2num(temp_STime(11:12))*60+...
           str2num(temp_STime(13:14));
       temp_ETime=datestr(sep_data(end).endTime,'yyyymmddHHMMss');
%        ETime=86400-(str2num(temp_ETime(9:10))*3600+str2num(temp_ETime(11:12))*60+...
%            str2num(temp_ETime(13:14)));
       Spoint=STime*sep_data(1).sampleRate;
%        Epoint=ETime*sep_data(end).sampleRate;
       if Spoint>0
           temp_merge.data=[zeros(Spoint,1);sep_data(1).data];% add zeros in the front
       end
       
      for k=1:length(sep_data)-1
          Mpoint=[];
          temp_T1=datestr(sep_data(k).endTime,'yyyymmddHHMMss');
          T1time=(str2num(temp_T1(9:10))*3600+str2num(temp_T1(11:12))*60+...
           str2num(temp_T1(13:14)));
          temp_T2=datestr(sep_data(k+1).startTime,'yyyymmddHHMMss');
          T2time=(str2num(temp_T2(9:10))*3600+str2num(temp_T2(11:12))*60+...
           str2num(temp_T2(13:14)));
          temp_T3=datestr(sep_data(k+1).endTime,'yyyymmddHHMMss');
          T3time=(str2num(temp_T3(9:10))*3600+str2num(temp_T3(11:12))*60+...
           str2num(temp_T3(13:14)));
          TTM=[T2time*sep_data(1).sampleRate:T3time*sep_data(1).sampleRate];
          
          Mpoint=zeros(round((T2time-T1time)*sep_data(1).sampleRate),1);
          temp_merge.data=[temp_merge.data;Mpoint];
          temp_merge.data=[temp_merge.data;sep_data(k+1).data];
      end
      Epoint=86400*sep_data(1).sampleRate-length(temp_merge.data);
      if Epoint>0
          temp_merge.data=[temp_merge.data;zeros(Epoint,1)];% add zeros in the end
      end
          temp_merge.startTime=datenum([temp_STime(1:8),'000000'],'yyyymmddHHMMss');
          temp_merge.endTime=datenum([temp_STime(1:8),'235959'],'yyyymmddHHMMss');
          Normallength=86400*sep_data(1).sampleRate;
          temp_merge.data=temp_merge.data(1:Normallength);
          temp_merge.sampleCount=Normallength;
    end


