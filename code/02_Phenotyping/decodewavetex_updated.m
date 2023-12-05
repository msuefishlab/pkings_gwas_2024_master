function eodtext=decodewavetex_updated(wave_text)
%DECODEWAVETEX
%
%   DECODEWAVETEX(WAVE_TEXT) retrieves information from the string
%   WAVE_TEXT, which is "information text" usually obtained from standard
%   Hopkin's Lab EOD files.  The returned argument is a structure
%   containing the following fields:
%       date            time            specimen        species
%       location        temperature     comment
%   These fields are populated by strings taken from parts of WAVE_TEXT.
%
%   See also READEODFILE, WRITEEODFILE

%   Written by B. Scott Jackson based partly on a previous version of
%   'decodewavetex.m' written by Christine Kinnier and Carl Hopkins
%   $Revision 1.00$  $Date: 2005/01/27 15:37$

% Remove null character from end
if wave_text(end) == 0
    wave_text(end) = [];
end

% Set default strings
eodtext.date=' ';
eodtext.time=' ';
eodtext.specimen=' ';
eodtext.species=' ';
eodtext.location=' ';
eodtext.temperature='25.0'; %if the temperature is not specified, set the default temperature to this value
eodtext.comment=wave_text;  %default all of text into comment

% Parse the 'wave_text' string
wtl = length(wave_text);
sc_stopvector = findstr(wave_text, ';');  %find semicolons in wave_text
space_stopvector = findstr(wave_text, ' ');


 %start BSJ Code
 
    if length(sc_stopvector) >= 3 && length(sc_stopvector)< 8   % if there are fewer than 3 semicolons in the text string, return defaults, else:
        if ~isempty(space_stopvector)
            stop1 = space_stopvector(1);
        else
            stop1 = wtl+1; % will always be greater than stop2 
        end
        
        stop2 = sc_stopvector(1);
        stop3 = sc_stopvector(2);
        stop4 = sc_stopvector(3);
	
        if length(sc_stopvector) >=4;   
            stop5 = sc_stopvector(4);   % if there are 4 or more semicolons, there was a comment
            eodtext.comment = wave_text(stop5+1:wtl);
        else
            stop5=  wtl+1;
            eodtext.comment = ' ';
        end;      
     
        if (stop1 < stop2) 
            eodtext.time = wave_text(1:stop1-1);
            eodtext.date = wave_text(stop1+1:stop2-1);
        else
            eodtext.time = ' ';
            eodtext.date = wave_text(1:stop2-1);
        end

        eodtext.specimen = wave_text(stop2+1:stop3-1);
        eodtext.species = wave_text(stop3+1:stop4-1);
        eodtext.location = wave_text(stop4+1:stop5-1);
	
        period_locs = strfind(eodtext.comment, '.');
        if ~isempty(period_locs) & (period_locs(end) == length(eodtext.comment))
            period_locs(end) = [];
        end
        for k = period_locs
            if (str2num(eodtext.comment(k-2:k+1)) >= 10.0) & (str2num(eodtext.comment(k-2:k+1)) <= 43.0)
                eodtext.temperature = eodtext.comment(k-2:k+1);
                break;
            end
        end
        
        %disp('this should not appear if this file was created before 2009')
        
    else
        %disp('this is a new 2009 eod file!')
        
        split_string=strsplit(wave_text, ';');
         
         tmp=strsplit(split_string{1},'=');
         eodtext.time = tmp{2};
         tmp=strsplit(split_string{2},'=');
         eodtext.date = tmp{2};
         tmp=strsplit(split_string{3},'=');
         eodtext.specimen = tmp{2};
         tmp=strsplit(split_string{4},'=');
         eodtext.species = tmp{2};
         tmp=strsplit(split_string{5},'=');
         eodtext.location = tmp{2};
         tmp=strsplit(split_string{6},'=');
         eodtext.temperature= tmp{2};
         tmp=strsplit(split_string{7},'=');
         eodtext.gain = tmp{2};
         tmp=strsplit(split_string{8},'=');
         eodtext.coupling = tmp{2};
         tmp=strsplit(wave_text,'Comments =');
         eodtext.comment = tmp{2};
        
    end
end
