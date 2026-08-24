%% FP_calib_analyze.m

clear; clc; close all;

%% ============================== 1. LOAD SAVED SCAN ==============================

[fname, fpath] = uigetfile('calibscan*.mat', 'Select a calibration scan file');
% if isequal(fname, 0)
%     error('No file selected.');
% end
FILE_TO_ANALYZE = fullfile(fpath, fname);

loaded = load(FILE_TO_ANALYZE, 'lambda', 'spectra', 'scan');
lambda = loaded.lambda;
spectra = loaded.spectra;
scan = loaded.scan;

folder = fpath;

fprintf('Loaded: %s\n', FILE_TO_ANALYZE);
fprintf('spectra size: %s, lambda size: %s\n', ...
    mat2str(size(spectra)), mat2str(size(lambda)));

%% ============================== 2. DIP ANALYSIS ==============================
% Per pixel: find the spectral dip (minimum). Build a 2D map of the dip
% wavelength over the scan grid, plus one sample spectrum so you can
% sanity-check what the dip pick is landing on.

lambda_nm = lambda(:) * 1e9;

dipVal_map = nan(scan.ny, scan.nx);
dipWl_map  = nan(scan.ny, scan.nx);

for ky = 1:scan.ny
    for kx = 1:scan.nx
        pixSpec = squeeze(spectra(kx, ky, :));
        [dp, dpIdx] = min(pixSpec);

        dipVal_map(ky, kx) = dp;
        dipWl_map(ky, kx)  = lambda_nm(dpIdx);
    end
end

fprintf('\n--- Analysis summary ---\n');
fprintf('Dip wavelength: mean %.3f nm, std %.3f nm\n', ...
    mean(dipWl_map(:)), std(dipWl_map(:)));

% --- Combined figure: reflectivity spectrum + dip wavelength map ---
kxC = round((scan.nx + 1) / 2);
kyC = round((scan.ny + 1) / 2);
centerSpec = squeeze(spectra(kxC, kyC, :));
[dpC, dpIdxC] = min(centerSpec);

figure('Name', 'Dip analysis', 'Position', [100 100 1000 400]);

subplot(1, 2, 1);
plot(lambda_nm, centerSpec); hold on;
plot(lambda_nm(dpIdxC), dpC, 'bv', 'MarkerFaceColor', 'b');
xlabel('Wavelength [nm]'); ylabel('Reflectivity');
legend('Spectrum', 'Dip'); grid on;

subplot(1, 2, 2);
imagesc(scan.x_um, scan.y_um, dipWl_map);
set(gca, 'YDir', 'normal'); axis image; colorbar;
set(gca, 'CLim', [1032 1033]);
xlabel('x [um]'); ylabel('y [um]'); title('Dip wavelength [nm]');


%% ============================== 3. YOKO-STYLE CALIBRATION FIGURE ==============================
% Self-contained re-implementation (normalize -> smooth -> derivative ->
% steepest-edge optimal wavelength). Tune the two settings below to match
% your spectra.

smoothSpectraWin = 10;    % spectrum smoothing window, in samples
smoothDerivWin   = 20;    % derivative smoothing window, in samples

% --- flatten the (nx,ny,nLambda) cube into (nLambda x Npts), one column per pixel ---
M = zeros(numel(lambda_nm), scan.Npts);
gridIdx = zeros(scan.Npts, 2);   % [kx, ky] for each column
col = 0;
for ky = 1:scan.ny
    for kx = 1:scan.nx
        col = col + 1;
        M(:, col) = squeeze(spectra(kx, ky, :));
        gridIdx(col, :) = [kx, ky];
    end
end

% --- normalize each spectrum to [0, 1], guarding against flat spectra ---
spanPerCol = max(M, [], 1) - min(M, [], 1);
spanPerCol(spanPerCol == 0) = eps;
Mnorm = (M - min(M, [], 1)) ./ spanPerCol;

% --- smooth, then compute + smooth the derivative ---
Msmooth = smoothdata(Mnorm, 1, 'movmean', smoothSpectraWin);
dM = zeros(size(Msmooth));
for col = 1:scan.Npts
    dM(:, col) = gradient(Msmooth(:, col), lambda_nm);
end
dMsmooth = smoothdata(dM, 1, 'movmean', smoothDerivWin);

% --- per-pixel optimal wavelength: steepest rising edge of the derivative ---
optimalWl_nm   = nan(1, scan.Npts);
derivAtOptimal = nan(1, scan.Npts);

for col = 1:scan.Npts
    [dVal, dIdx] = max(dMsmooth(:, col));
    optimalWl_nm(col)   = lambda_nm(dIdx);
    derivAtOptimal(col) = dVal;
end

% --- reshape per-pixel results into 2D maps ---
optimalWlMap = nan(scan.ny, scan.nx);
derivMap     = nan(scan.ny, scan.nx);
for col = 1:scan.Npts
    kx = gridIdx(col, 1); ky = gridIdx(col, 2);
    optimalWlMap(ky, kx) = optimalWl_nm(col);
    derivMap(ky, kx)     = derivAtOptimal(col);
end

% --- Smoothing/filtering check: raw vs smoothed spectrum, center pixel ---
kxC = round((scan.nx + 1) / 2);
kyC = round((scan.ny + 1) / 2);
centerCol = find(gridIdx(:, 1) == kxC & gridIdx(:, 2) == kyC, 1);

figure('Name', 'Smoothing check');
plot(lambda_nm, Mnorm(:, centerCol), 'Color', [0.7 0.7 0.7]); hold on;
plot(lambda_nm, Msmooth(:, centerCol), 'b', 'LineWidth', 1.5);
xlabel('Wavelength [nm]'); ylabel('Normalized intensity');
title(sprintf('Center pixel (kx=%d, ky=%d): raw vs smoothed/filtered', kxC, kyC));
legend('Raw (normalized)', 'Smoothed/filtered'); grid on;

% --- Main figure: 2 panels ---
figure('Name', 'Calibration analysis (yoko style)', 'Position', [50 50 1000 450]);

subplot(1, 2, 1);
imagesc(scan.x_um, scan.y_um, derivMap);
set(gca, 'YDir', 'normal'); axis image; colormap(gca, 'hot'); colorbar;
xlabel('x [um]'); ylabel('y [um]'); title('Max derivative map');

subplot(1, 2, 2);
imagesc(scan.x_um, scan.y_um, optimalWlMap);
set(gca, 'YDir', 'normal'); axis image; colormap(gca, 'turbo'); colorbar;
set(gca, 'CLim', [1033.1 1033.5]);
xlabel('x [um]'); ylabel('y [um]'); title('Optimal wavelength map [nm]');

try
    saveas(gcf, fullfile(folder, 'calib_analysis_yoko_style.png'));
catch
end

fprintf('\nAnalysis complete.\n');