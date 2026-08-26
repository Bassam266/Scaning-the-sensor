%% FP_calib_dipmap_analyze.m
% Load a saved FP calibration scan, smooth the spectra (optionally over a
% restricted wavelength range) using the FPISpectraAnalysis class, then
% plot the processed reflectivity spectra together with the resulting
% dip-position map.
%
%   1. Load the saved scan file
%   2. Crop to [LAMBDA_MIN_NM, LAMBDA_MAX_NM] (optional) and smooth over
%      SMOOTH_WIN points, using FPISpectraAnalysis
%   3. Plot the processed reflectivity spectra and the dip-position map
%
% Requires on path: FPISpectraAnalysis.m

clear; clc; close all;

%% ============================== 1. LOAD SAVED SCAN ==============================

[fname, fpath] = uigetfile('calibscan*.mat', 'Select a calibration scan file');
if isequal(fname, 0)
    error('No file selected.');
end
FILE_TO_ANALYZE = fullfile(fpath, fname);

loaded = load(FILE_TO_ANALYZE, 'lambda', 'spectra', 'scan');
lambda = loaded.lambda;
spectra = loaded.spectra;
scan = loaded.scan;

fprintf('Loaded: %s\n', FILE_TO_ANALYZE);
fprintf('spectra size: %s, lambda size: %s\n', ...
    mat2str(size(spectra)), mat2str(size(lambda)));

%% ============================== 2. SETTINGS ==============================

SMOOTH_WIN = 15;   % moving-average window, in samples, for Smooth()

% Wavelength range to keep, in nm. Leave either bound as [] to use the
% full range on that side, e.g. LAMBDA_MIN_NM = 1032; LAMBDA_MAX_NM = [];
LAMBDA_MIN_NM = [];
LAMBDA_MAX_NM = [];

% Colorbar limits for the dip-position map. Leave as [] to auto-scale
% to the map's own min/max, or set explicit [min max] values, e.g. [1032 1034].
MAP_CLIM = [];

%% ============================== 3. CROP TO WAVELENGTH RANGE ==============================

lambda_nm = lambda(:) * 1e9;

keepMask = true(size(lambda_nm));
if ~isempty(LAMBDA_MIN_NM)
    keepMask = keepMask & (lambda_nm >= LAMBDA_MIN_NM);
end
if ~isempty(LAMBDA_MAX_NM)
    keepMask = keepMask & (lambda_nm <= LAMBDA_MAX_NM);
end
if ~any(keepMask)
    error('No wavelength points left after applying LAMBDA_MIN_NM/LAMBDA_MAX_NM.');
end

lambdaCropped = lambda(keepMask);
spectraCropped = spectra(:, :, keepMask);

fprintf('Using %d/%d wavelength points (%.3f-%.3f nm) after cropping.\n', ...
    nnz(keepMask), numel(keepMask), min(lambda_nm(keepMask)), max(lambda_nm(keepMask)));

%% ============================== 4. LOAD INTO FPISpectraAnalysis & SMOOTH ==============================
% spectraCropped is (nx, ny, nLambda). The class's constructor only
% avoids a lambdaDim bug when a *trailing* size-1 dimension disappears
% during permute, so add back a dummy dimension before nLambda and use
% lambdadim=4 (see FP_calib_peaks_analyze.m for the same workaround).
spectra4D = reshape(spectraCropped, [scan.nx, scan.ny, 1, numel(lambdaCropped)]);

fpi = FPISpectraAnalysis('lambda', lambdaCropped, 'spectra', spectra4D, 'lambdadim', 4);

fpi.Smooth(SMOOTH_WIN, 'spectra');

%% ============================== 5. DIP-POSITION MAP ==============================
% fpi.spectra is (nLambda, nx, ny) after the constructor's permute.

[dipVal_xy, dipIdx_xy] = min(fpi.spectra, [], 1);
dipVal_xy = squeeze(dipVal_xy);              % nx x ny
dipIdx_xy = squeeze(dipIdx_xy);              % nx x ny
dipWl_xy_nm = fpi.lambda(dipIdx_xy) * 1e9;   % nx x ny

dipVal_map = dipVal_xy.';       % ny x nx, matches imagesc(scan.x_um, scan.y_um, ...)
dipWl_map  = dipWl_xy_nm.';     % ny x nx

fprintf('\n--- Dip position summary (after smoothing) ---\n');
fprintf('Dip wavelength: mean %.3f nm, std %.3f nm\n', ...
    mean(dipWl_map(:), 'omitnan'), std(dipWl_map(:), 0, 'omitnan'));

%% ============================== 6. PLOTS ==============================
% Both panels together in a single figure.

figure('Name', 'Spectra & dip map', 'Position', [100 100 1000 450]);

% --- left: processed reflectivity spectra, all measurement points overlaid ---
subplot(1, 2, 1);
plot(fpi.lambda * 1e9, fpi.spectra(:, :));
xlabel('Wavelength [nm]'); ylabel('Reflectivity');
title('Reflectivity spectra (smoothed)');
grid on; box on;

% --- right: dip-position map ---
subplot(1, 2, 2);
imagesc(scan.x_um, scan.y_um, dipWl_map);
set(gca, 'YDir', 'normal'); axis image; colormap(gca, 'turbo');
cb = colorbar;
cb.Label.String = 'Dip wavelength [nm]';
if ~isempty(MAP_CLIM)
    set(gca, 'CLim', MAP_CLIM);
end
xlabel('x [um]'); ylabel('y [um]');
title('Dip position map (after smoothing)');

fprintf('\nDone.\n');
