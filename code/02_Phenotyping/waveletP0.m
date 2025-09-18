%% ---------- CONFIG ----------
N_target   = 256;     % resample each segment to this many points
maxK       = 8;       % explore k = 2..maxK for k-means
use_cwt    = true;    % set false to use DWT fallback
waveletFun = 'amor';  % CWT mother ('amor' = analytic Morse)

%% ---------- 1) Extract segments (iZC1 - width ... iZC1) ----------
n = numel(averaged_eods);
segments = cell(n,1);

for i = 1:n
         measurement_data(i)=standard_eod_measurement(averaged_eods(i).wave,averaged_eods(i).sampRate,averaged_eods(i).sample_name);

    sr   = averaged_eods(i).sampRate;
    x    = averaged_eods(i).wave(:);                 % ensure column
    width = fix(sr * 0.0005);                        % your 0.5 ms window
    i1   = measurement_data(i).iZC1 - width;
    i2   = measurement_data(i).iZC1;

    % bounds check
    i1 = max(1, i1);
    i2 = min(numel(x), i2);

    seg = x(i1:i2);

    % If the segment ended up too short (edge cases), pad symmetrically with endpoints
    if numel(seg) < 2
        seg = repmat(seg, max(2, numel(seg)), 1);
    end

    % Normalize amplitude (z-score) to focus on shape, not scale
    seg = seg - mean(seg);
    ssd = std(seg);
    if ssd > 0
        seg = seg / ssd;
    end

    % Resample to common length
    segments{i} = resample(seg, N_target, numel(seg));
end

% Stack to matrix (samples x observations)
X = cell2mat(segments');     % size = N_target x n
X = X';                      % now n x N_target (rows = observations)

%% ---------- 2) Wavelet features ----------
if use_cwt
    % Continuous Wavelet Transform; summarize power across time → scale spectrum
    % This gives translation-robust spectral-shape features.
    % Wavelet Toolbox required for cwt.
    all_feats = [];
    for i = 1:n
        sig = X(i, :).';
        % CWT returns coefficients as (numScales x numTime)
        [cfs, freqs] = cwt(sig, waveletFun, N_target);   %#ok<ASGLU>
        pow = abs(cfs).^2;
        scale_spectrum = mean(pow, 2);                   % average over time
        % log-compress
        scale_spectrum = log10(scale_spectrum + eps);
        all_feats(:, i) = scale_spectrum;                %#ok<SAGROW>
    end
    F = all_feats.';   % n x numScales
else
    % ---- Fallback using discrete wavelet packet-ish summary via DWT ----
    % No Wavelet Toolbox? You can still use wavedec (requires Wavelet Toolbox) or
    % build your own filterbank. Below uses wavedec if available.
    % (If wavedec isn’t available, replace this section with FFT-band powers.)
    lvl = 4; wname = 'db4';
    % pre-allocate using first signal to get coeff length
    [c0, l0] = wavedec(X(1,:).', lvl, wname);
    F = zeros(n, numel(c0));
    F(1,:) = c0(:).';
    for i = 2:n
        [ci, ~] = wavedec(X(i,:).', lvl, wname);
        F(i,:) = ci(:).';
    end
    % Standardize DWT coeffs
    F = (F - mean(F,1)) ./ (std(F,[],1) + 1e-8);
end

%% ---------- 3) Dimensionality reduction (PCA) ----------
[Fcoeff, Fscore, ~, ~, explained] = pca(F, 'Centered', true);

% Choose number of PCs to cover ~95% variance (or cap at 10 for plotting convenience)
cumexp = cumsum(explained);
numPCs = find(cumexp >= 95, 1, 'first');
numPCs = min(max(numPCs, 2), 10);
Z = Fscore(:, 1:numPCs);   % embedding for clustering

%% ---------- 4) Pick k via silhouette and cluster ----------
bestK = 2;
bestSil = -Inf;
sil_cache = cell(0);

opts = statset('UseParallel', false, 'MaxIter', 1000, 'Display', 'off');

for k = 2:maxK
    % multiple restarts for robustness
    [idx_k, ~] = kmeans(Z, k, 'Replicates', 10, 'MaxIter', 1000, ...
                        'Options', opts, 'Distance', 'sqeuclidean');
    s = silhouette(Z, idx_k);
    meanSil = mean(s);
    sil_cache{end+1} = struct('k',k,'idx',idx_k,'sil',s,'meanSil',meanSil); %#ok<SAGROW>
    if meanSil > bestSil
        bestSil = meanSil;
        bestK = k;
    end
end

% Refit final k-means at chosen k (or use the cached one with best silhouette)
chosen = sil_cache{cellfun(@(a) a.k==bestK, sil_cache)};
idx = chosen.idx;

fprintf('[Clustering] Best k = %d (mean silhouette = %.3f)\n', bestK, chosen.meanSil);

%% ---------- 5) Visualizations ----------
% 2D PCA scatter
figure('Name','PCA scatter (first 2 PCs)');
gscatter(Fscore(:,1), Fscore(:,2), idx);
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title(sprintf('CWT+PCA embedding colored by k-means clusters (k=%d)', bestK));
grid on; box on;

% Silhouette plot
figure('Name','Silhouette (chosen k)');
silhouette(Z, idx);
title(sprintf('Silhouette plot for k = %d (mean = %.3f)', bestK, chosen.meanSil));

% Cluster-average waveforms (in the resampled segment domain)
figure('Name','Cluster-average segments');
window_sec = 0.0005;            % duration of your segment in seconds
t = linspace(-window_sec, 0, N_target);
for k = 1:bestK
    subplot(bestK,1,k);
    members = (idx == k);
    if any(members)
        meanSeg = mean(X(members,:), 1);
        plot(t*1e3, meanSeg, 'LineWidth', 1.5); hold on;
        plot(t*1e3, X(members,:).', 'Color', [0.7 0.7 0.7 0.3]); % faint members
        hold off;
        ylabel(sprintf('k = %d', k));
    else
        text(0.5, 0.5, sprintf('No members in k=%d',k), 'HorizontalAlignment','center');
    end
    if k == 1
        title('Cluster means (bold) with member segments (faint)');
    end
    if k == bestK
        xlabel('Time relative to iZC1 (ms)');
    end
    grid on; box on;
end

%% ---------- 6) Optional: label & inspect exemplars ----------
% Find closest-to-centroid exemplar per cluster
exemplars = zeros(bestK,1);
for k = 1:bestK
    members = find(idx == k);
    if isempty(members), continue; end
    c_k = mean(Z(members,:), 1);
    [~, j] = min(sum((Z(members,:) - c_k).^2, 2));
    exemplars(k) = members(j);
end

disp('Exemplar indices per cluster (rows):');
disp(exemplars(:).');

% Example: view the scalogram of an exemplar (requires Wavelet Toolbox)
if use_cwt && any(exemplars)
    k1 = find(exemplars, 1, 'first');
    i_ex = exemplars(k1);
    figure('Name','Scalogram of an exemplar');
    cwt(X(i_ex,:).', waveletFun, N_target);
    title(sprintf('Cluster %d exemplar (index %d)', k1, i_ex));
end

%% ---------- Notes ----------
% • CWT features: we used the scale-spectrum (mean power across time) to make
%   features translation-robust. You can also concatenate time-marginals or
%   use more scales via 'VoicesPerOctave' in cwt for finer resolution.
% • If segments vary in sampling rate, our resampling step standardizes them.
% • If you prefer not to decide k, try hierarchical clustering:
%     Y = pdist(Z); ZL = linkage(Y,'ward'); figure; dendrogram(ZL);
%   then “cut” the tree at different heights (cluster(ZL, 'Maxclust', k)).
% • Outlier-robust option: try DBSCAN on Z:
%     eps = 0.8; MinPts = 5; labels = dbscan(Z, eps, MinPts);
%   points labeled -1 are noise/outliers.
