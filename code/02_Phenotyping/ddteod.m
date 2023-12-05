function [ddtwave]=ddteod(wave,s_rate)
%%%%%%%%%%
%ddteod.m
%first derivative of wave is estimated using a numerical method and plotted
%%%%%%%%%%
%
%call this function after normalizing wave with function normeodp1.m
%vector 'wave' (n_pts x 1) scaled voltage data (sitting on 0 baseline and with peak-to-peak of 1)
%vector 'time' (1 x n_pts) centered time data (in msec, such that tP1=0 is the reference)
%
%this function works on a filtered wave (to reduce high frequency noise to enable it to find max slopes
%in minute waveform phases
%
%numerical estimation of first derivative sums triangle-weighted, pairwise differences of specified #
%(smooth/2) of points before index and same # (smooth/2) points after index (method used by maf in
%program 'automeas')
%
[~,n_pts]=size(wave);
%--set parameters & initialize---------
smooth = 8;										%how many pnts to sample for calc ddt? (performs best between 6 & 16)
ddtwave = zeros((size(wave)));			%initialize new vect (n_pts x 1) to store deriv of wave
%--filter out high freq noise----------%same code as filteod.m
freq1=1;											%low cutoff freq in Hz
freq2=s_rate/8;								%high cuttoff freq in Hz (1/4 x the Nyquist or 1/8 x sampling freq)
fs=[freq1 freq2]*2/s_rate;
[b,a]=butter(3,fs);							%third order butterworth digital filter
filtwave=filter(b,a,wave);
%--approx 1stDeriv:triangle weighting--
factor = (4/(smooth*(smooth+2)+(rem(smooth,2))));
for i = 1:+1:(n_pts-smooth)
   for j = 0:+1:(smooth-1)
      if j<(smooth/2)
         k=j+1;
         else k=smooth-j;
         end
    	ddtwave((i+(smooth/2)))=cumsum(((filtwave(i+j+1)-filtwave(i+j))*factor*k));
      end
   end
for i = 1:+1:(round(smooth/2))			%clean up edge effects: pad the first & last parts of the derivative
   ddtwave(i) = ddtwave(((smooth/2)+1));
	end
for i = 1:+1:((smooth/2)-1+(rem(smooth,2)))
   ddtwave((n_pts-i))=ddtwave((n_pts-(smooth/2)));
   end
t1 = -3;
t2 = 3;
ddtwave = s_rate*ddtwave;			%scaling needed so all slopes on same time base