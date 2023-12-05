function [vT1,tT1,iT1,vT2,tT2,iT2]=eod_start_end(wave,time,thresh)
%%%%%%%%%%
%Find EOD Start and End
%Finds the beginning and end of an EOD waveform, defined by waveform crossing
%threshold 'thresh'.  Must specify a sample rate.  Returns
% vT1 - Voltage at Start
% tT1 - Time of Start
% iT1 - Index of Start
% vT2 - Voltage at End
% tT2 - Time of End
% iT2 - Index of End

[tmin,istart] = min(time);			%first index of eod
[tmax,iend] = max(time);			%last index of eod

for j=istart:+1:iend			%look for first point that deviates >= 2% ptp
   if abs(wave(j))>thresh
            if abs(wave(j+1))>thresh
                iT1 = j;								%index of T1
                tT1 = time(iT1);						%time of T1
                vT1 = wave(iT1);						%voltage at T1
                break
            else
                j=j+1;
                continue
            end
    else
        j=j+1;
        continue
   end
end

k = iend;
while abs(wave(k))<thresh			%read backwards to find T2 (2% deviation)
   k = k - 1;
end
iT2 = k;									%index of T2
tT2 = time(k);							%time of T2
vT2 = wave(k);		