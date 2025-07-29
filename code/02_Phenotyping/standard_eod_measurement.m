function [measurement_data] = standard_eod_measurement(wave,s_rate,samp_name)
%%%%%%%%%%
%STANDARD_EOD_MEASUREMENT
%performs EOD measurements on previously normalized EOD waveform, based on
%measurescript.m
%modified 3 November 2008, for JRG's study of signal variation in P.
%kingsleyae, based on MEA's code for analyzing variation in P. magnostipes
%reproductive isolation.

%%--Settings: Should be Paramaterized...---
threshold = 0.02;
firstforty = wave(1:40);
stddev = std(firstforty);
if stddev > threshold				%if wave too noisy
   threshold = 2*stddev;			%set noise threshold = 2 x std dev
end
[~,n_pts]=size(wave);
time=linspace(0,1000*n_pts/s_rate,n_pts);

%% --calc first derivative used for calculating slopes---------
ddtwave=ddteod(wave,s_rate);									%returns vector 'ddtwave'

%% --Describe Peak P1, always present, by definition ------ %%
[~,measurement_data.iP1]=max(wave);
time=time-time(measurement_data.iP1);
measurement_data.tP1 = time(measurement_data.iP1);                                          %time of P1
measurement_data.vP1 = wave(measurement_data.iP1);                                          %voltage at P1
[measurement_data.sS1,measurement_data.iS1] = max(ddtwave(1:measurement_data.iP1));         %index (iS1) and Slope (sS1) of Maximum Slope of P1 (S1)
measurement_data.tS1 = time(measurement_data.iS1);                                          %Time of S1
measurement_data.vS1 = wave(measurement_data.iS1);                                          %voltage at S1

%% --Describe Beginning and End of Waveform, always present ---- %%
[measurement_data.vT1,measurement_data.tT1,measurement_data.iT1,measurement_data.vT2,measurement_data.tT2,measurement_data.iT2]=eod_start_end(wave,time,threshold);

%% --Describe P0.  Always Assumed to be present -------- %%
for k = measurement_data.iP1:-1:1				
	if wave(k) < 0	       
        measurement_data.iZC1 = k;				
		break;
    else
        measurement_data.iZC1 = measurement_data.iP1;
    end
end

measurement_data.tZC1 = time(measurement_data.iZC1);				%time of ZC1
measurement_data.vZC1 = wave(measurement_data.iZC1);				%voltage at ZC1
measurement_data.sZC1 = ddtwave(measurement_data.iZC1);             %slope at ZC1 

width=fix(s_rate*.0005);
markpoint=measurement_data.iZC1-width;
measurement_data.aP0=sum(wave(markpoint:measurement_data.iZC1))/s_rate;
measurement_data.aP0=measurement_data.aP0/1e-6;

[measurement_data.vP0,measurement_data.iP0] = min(wave(markpoint:measurement_data.iZC1)); %attempt to find the minimum
measurement_data.tP0 = time(measurement_data.iP0+markpoint);     %attempt to find the time of minimum

%% --Describe Peak P2, assumed always present ------ %%
[measurement_data.vP2,measurement_data.iP2]=min(wave(measurement_data.iP1:n_pts));                 %Voltage and Index of P2
measurement_data.iP2 = measurement_data.iP1 + measurement_data.iP2 - 1;                                 %Adjust index of P2
measurement_data.tP2 = time(measurement_data.iP2);                                                                       %time of P2
[measurement_data.sS2,measurement_data.iS2] = min(ddtwave(measurement_data.iP1:measurement_data.iP2));  %index (iS2) and Slope (sS2) of Maximum Slope of P1-P2 (S2)
measurement_data.iS2 = measurement_data.iP1 + measurement_data.iS2 - 1;                                                  %Adjust index of S2
measurement_data.tS2 = time(measurement_data.iS2);                                                                                        %Time of S2
measurement_data.vS2 = wave(measurement_data.iS2);                                                                                        %Voltage of S2

%% --Zero Crossings ------ %%
dumpos1 = n_pts-1;						%default near end of wave
dumpos2 = n_pts;						%default at end of wave
for i = measurement_data.iP1:+1:n_pts					%start at P1 and read forward
	if wave(i) <= threshold		
   	dumpos1 = i - 1;				%index just before wave crosses +threshold
		break;     						%break out of for loop and leave dumpos1 = corrected i
  	end
end
for i = dumpos1:+1:n_pts				%start at dumpos1 and read forward
	if wave(i) <= (-threshold)	
     	dumpos2 = i;					%index where wave crosses -threshold
     	dumpos2_present = 1;			%if so, set cond for phase2 presence to true      
     	break;    						%break out of for loop and leave dumpos2 = new i
    end
end

if dumpos2_present					%if condition true, find ZC2 landmarks
   j = round((dumpos1 + dumpos2)/2);
   if abs(wave(dumpos1))<abs(wave(j))
      measurement_data.iZC2=dumpos1;
   elseif abs(wave(dumpos2))<abs(wave(j))
     	measurement_data.iZC2=dumpos2;					%these statements check for rounding error
   else									%and find the actual voltage nearest zero-xing
      measurement_data.iZC2=j;							%in this (usually) steep voltage transition
   end
  	measurement_data.tZC2 = time(measurement_data.iZC2);				%time of ZC2
  	measurement_data.vZC2 = wave(measurement_data.iZC2);				%voltage at ZC2
  	measurement_data.sZC2 = ddtwave(measurement_data.iZC2);			%slope at ZC2
end

%% --calculate phase areas---------
measurement_data.aP1 = (sum(wave(measurement_data.iT1:measurement_data.iZC2)))/s_rate;
measurement_data.aP2 = (sum(wave(measurement_data.iZC2:measurement_data.iT2)))/s_rate;

%% --calculate the power spectrum-------
% PS=power_spectrum_eod(wave,s_rate);	%call to func that plots power spectr & calcs fftstats
% 
% measurement_data.Fmax = PS.Fmax;
% measurement_data.Pmax = PS.Pmax;
% measurement_data.Flow = PS.Flow;
% measurement_data.Plow = PS.Plow;
% measurement_data.Fhigh = PS.Fhigh;
% measurement_data.Phigh = PS.Phigh;
% measurement_data.Bndwidth = PS.Fhigh - PS.Flow;
measurement_data.sample_name=string(samp_name);
