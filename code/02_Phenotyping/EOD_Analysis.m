root=get_root();

tsvFilename = fullfile(root,'input_data/09_RNASeq/rnaseq_metadata.txt'); % Replace with your actual TSV file name
[subjectlist,samplenamelist] = getSubjectListFromTSV(tsvFilename);

normalized_eods = {};
averaged_eods = {};
a = 1;

for b = 1:length(subjectlist)
    if subjectlist{b} == "ND"
        continue
    elseif isempty(subjectlist{b})
        continue
    else
        filename = fullfile(root,'input_data/09_RNASeq/rnaseq_eods',subjectlist{b}); % Assuming full path is provided in the TSV
        samplename=samplenamelist{b};
        eod = loadEODData(filename);
        if isempty(eod)
            continue;
        end
        
        message = ['Currently working on ', samplename, ' from file ', filename];
        disp(message);
        neodwave = processEODWaveforms(eod);
        normalized_eods = updateEODStruct(normalized_eods, neodwave, eod, samplename, a);
        averaged_eods = updateAveragedEODStruct(averaged_eods, neodwave, eod, samplename, a);
        a = a + 1;
    end
end

function gitRoot=get_root()
    [status, cmdout] = system('git rev-parse --show-toplevel');
    if status == 0
        gitRoot = strtrim(cmdout);
    else
        disp('Not a Git repository or Git is not installed.');
    end
end

function [subjectlist,samplenamelist] = getSubjectListFromTSV(tsvFilename)
    opts = detectImportOptions(tsvFilename, 'FileType', 'text');
    tsvData = readtable(tsvFilename, opts);
    subjectlist = tsvData.EOD_FILE;
    samplenamelist=tsvData.PN;
end


function eod = loadEODData(filename)
    [~, ~, extension] = fileparts(filename);
    if strcmpi(extension, '.eod')
        eod = ReadEODFileWithConversionto2019(filename);
    else
        load(filename);
        eod = updateEODFormat(eod);
    end
end

function eod = updateEODFormat(eod)
    % Try to convert the date using a standard format
    try
        eodDate = datetime(datestr(eod(1).date));
    catch
        % If the standard format fails, try the 'July 06, 2019' format
        try
            eodDate = datetime(eod(1).date, 'InputFormat', 'MMMM dd, yyyy');
            disp(eodDate)
        catch
            % Handle other potential errors or formats
            disp('Unknown date format');
            return;
        end
    end

    if eodDate < datetime('1/1/2019','InputFormat','MM/dd/uuuu')
        % Old format EOD, convert fields to proper names for compatibility
        [eod.('specimenno')] = eod.('specimen');
        [eod.('temp')] = eod.('temperature');
        [eod.('Rate')]=eod.('s_rate');
        [eod.('Range')]=eod.('adrange');
    else
        % New format, do nothing
    end
end


function neodwave = processEODWaveforms(eod)
        neodwave={};
        % for each waveform in the EOD file
        for i=1:size(eod,2)
            %initialize variables
            pospeakidx=[];
            negpeakidx=[];
            pospeakidx_new=[];
            negpeakidx_new=[];

            %calculate baseline
            baseline = mean(eod(i).wave(1:100));
            
            % filter out very long or very short recordings (spurrious)
            if size(eod(i).wave,1)> 5000 || size(eod(i).wave,1) < 1000
                continue
            else

                %find the maximum and minimum values
                [maxvalue, pospeakidx]=max(eod(i).wave);
                [minvalue,negpeakidx]=min(eod(i).wave);
                
                %find the indicies of those maxima
                IMAX=find(eod(i).wave==maxvalue);
                IMIN=find(eod(i).wave==minvalue);
                
                % filter out Clipped EODs
                if size(IMAX,1)> 2 || size(IMIN,1)> 2
                    disp("Clipped EOD detected!")
                    continue
                else
                    %standardize orientation so that the first peak in the
                    %recording is positive
                    if pospeakidx > negpeakidx 
                        eod(i).wave=-eod(i).wave;
                        [maxvalue, pospeakidx]=max(eod(i).wave);
                        [minvalue,negpeakidx]=min(eod(i).wave);
                    end
                end
                    
                %normalize p-t-p amplitude
                    neodwave{i} = eod(i).wave-baseline;
                    neodwave{i} = neodwave{i}/(neodwave{i}(pospeakidx)-neodwave{i}(negpeakidx));

                % center waveforms so that P1 occurs at t=0
                 offset=(size(neodwave{i},1)/2)-pospeakidx  ;
                 neodwave{i}=circshift(neodwave{i},floor(offset),1);
                
                fprintf('Processing waveform %d: Old Baseline = %f New Baseline= %f\n', i, baseline,mean(neodwave{i}(1:100)));

            end
        end

end

function normalized_eods = updateEODStruct(normalized_eods, neodwave, eod, sample_name, a)

      neodwave_filt=filterByMostCommonDimension(neodwave);

      normalized_eods(a).sample_name=sample_name;
      normalized_eods(a).version=1;
      normalized_eods(a).date=eod.date;
      normalized_eods(a).time=eod.time;
      normalized_eods(a).specimen=eod.specimenno;
      normalized_eods(a).species=eod.species;
      normalized_eods(a).location=eod.location;
      normalized_eods(a).temperature=eod.temp;
      normalized_eods(a).comment=eod.comments;
      normalized_eods(a).sampRate=eod.Rate;
      normalized_eods(a).ADrange=eod.Range;
      normalized_eods(a).wave = cell2mat(neodwave_filt);
end

function averaged_eods = updateAveragedEODStruct(averaged_eods, neodwave, eod, sample_name, a)
      neodwave_filt=filterByMostCommonDimension(neodwave);

      M=cat(2,neodwave_filt{:});
      averaged_eods(a).sample_name=sample_name;
      averaged_eods(a).version=1;
      averaged_eods(a).date=eod.date;
      averaged_eods(a).time=eod.time;
      averaged_eods(a).specimen=eod.specimenno;
      averaged_eods(a).species=eod.species;
      averaged_eods(a).location=eod.location;
      averaged_eods(a).temperature=eod.temp;
      averaged_eods(a).comment=eod.comments;
      averaged_eods(a).sampRate=eod.Rate;
      averaged_eods(a).ADrange=eod.Range;
      
      [~,h]=size(M);
      if h > 1
         averaged_eods(a).wave=mean(M');
      else
          averaged_eods(a).wave=M';
      end
end

function filteredArray = filterByMostCommonDimension(cellArray)
    % Initialize a map to store the count of each dimension
    dimCount = containers.Map('KeyType', 'char', 'ValueType', 'int32');

    % Iterate through the cell array and count dimensions
    for i = 1:length(cellArray)
        dims = size(cellArray{i});
        dimKey = sprintf('%dx%d', dims(1), dims(2)); % Convert dimensions to a string key
        if isKey(dimCount, dimKey)
            dimCount(dimKey) = dimCount(dimKey) + 1;
        else
            dimCount(dimKey) = 1;
        end
    end

    % Find the most common dimension, preferring longer dimensions in case of a tie
    maxCount = 0;
    mostCommonDim = '';
    dimKeys = keys(dimCount);
    for i = 1:length(dimKeys)
        currentCount = dimCount(dimKeys{i});
        if currentCount > maxCount || ...
           (currentCount == maxCount && getLongerDimension(dimKeys{i}, mostCommonDim))
            maxCount = currentCount;
            mostCommonDim = dimKeys{i};
        end
    end

    % Filter the cell array to include only elements with the most common dimension
    filteredArray = {};
    for i = 1:length(cellArray)
        currentDim = sprintf('%dx%d', size(cellArray{i}, 1), size(cellArray{i}, 2));
        if strcmp(currentDim, mostCommonDim)
            filteredArray{end+1} = cellArray{i};
        end
    end
end

function isLonger = getLongerDimension(dim1, dim2)
    % Helper function to determine if the first dimension is longer than the second
    dims1 = sscanf(dim1, '%dx%d');
    dims2 = sscanf(dim2, '%dx%d');
    isLonger = prod(dims1) > prod(dims2);
end


