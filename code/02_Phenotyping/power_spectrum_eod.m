function [PS] = power_spectrum_eod(x,Fs)
%FFTWAVE Power Spectrum of EOD Waveform
%
%	Syntax [PS,FFTSTATS] = FFTWAVE(x,Fs) 
%
%	Uses the matlab function pwelch , which
%	estimates the power spectrum of a discrete-time signal vector
%	using Welch's averaged modified periodogram method.  A Hanning
%	window is employed .  Calculates the power spectrum in
%	terms of dB relative to max power vs log10(freq) in figure(3).
%
%	Takes the input variables:
%		x (digitized waveform: 'wave') 
%		Fs (sampling frequency in Hz: 's_rate')
%
%	Returns the 1*3 vector fftstats, the elements of which are:
%		(1,1) Flow (lower frequency 3dB down from max)
%		(1,2) Fmax (frequency with max power)
%		(1,3) Fhigh (higher frequency 3dB down from max)
%PW


% see https://www.mathworks.com/matlabcentral/answers/18617-pwelch-vs-psd
% Hi Kevin,
% 
% The result of psd is not correctly scaled, therefore you should use pwelch. For spectral density, the result should be scaled by the sampling frequency, which is not performed by psd.
% 
% If you look at the two results, the f vector should be the same. If you take the ratio of two pxx, you can that most samples, the differ by a factor of 125, which is basically half of sampling frequency. This is because the resulting spectrum is one-sided. The two end points differs by a factor of 250. This is due to the fact that these two points correspond to DC and Nyquist frequency and should not be doubled even if it's one-sided.
% 
% BTW, for your case, there is really no overlap happening because your data and window length are the same.
% 
% HTH
% Thank you very much Honglei for your response. However, I have a question about your last comment about there being no overlap. You said the data and the window length are the same, but I never explicitly said how long the data was. Did you mean that there is no overlap because the nfft is the same size as the window? Or, did you think that my variable "data" is the same length my window? The length of "data" is actually at least 10,000 samples long. I just want to make sure that I am using the correct parameter values to make sure that my samples have 50% overlap.
% 
% Thanks again!
% N = 1024;
% 
% Here's a small tweak to Rick's code which takes into consideration the overlap value and the fact that the DC and Nyquist values should not be scaled. Also, I scaled the results of psd.m instead of pwelch.m. Now Pxx and Pyy are identical.
% 
%  load handel;
%  N = 1024;   
%  [Pxx,Fx] = psd(y,N,Fs,hanning(N));   
%  Pxx(2:end-1) = Pxx(2:end-1)*2;
%  Pxx = Pxx/Fs;  
%  [Pyy,Fy] = pwelch(y, hanning(N),0, N, Fs);
%  plot(Fx,10*log10(Pxx),Fy,10*log10(Pyy));
%  legend('psd','pwelch');
% -Paul

%[P,F] = psd(x,length(x),Fs); %deprecated
[PS.P,PS.F]=pwelch(x, hanning(length(x)),0, length(x), Fs);
%[pxx,f] =   pwelch(x, window,          noverlap, f,  fs) 

PS.P = PS.P/max(PS.P);
PS.PdB = 10*log10(PS.P);						%scale the power spectrum in dB relative to max
[Maxp,imax] = max(PS.P);					%find the index of max power
for j = (imax:-1:1)						%read from frequency w/ max power toward lower frequencies
   if PS.PdB(j) <= -3						%find the frequency at which power is 3 dB down from max
   	ilow = j;
   	break;
   else
      ilow = 1;							%if it never reaches -3dB, set low frequency to DC
   end
end
for k = (imax:+1:size(PS.F))				%read from frequency w/ max power toward higher frequencies
   if PS.PdB(k) <= -3						%find the frequency at which power is 3 dB down from max
   	ihigh = k;							%will always find, so else not necessary
   	break;
   end
end
PS.Fmax = PS.F(imax);
PS.Pmax = PS.PdB(imax);
PS.Flow = PS.F(ilow);
PS.Plow = PS.PdB(ilow);
PS.Fhigh = PS.F(ihigh);
PS.Phigh = PS.PdB(ihigh);
PS.Bndwidth = PS.Fhigh - PS.Flow;