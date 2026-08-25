%% FP_calib_analyze.m
% Analysis of a saved FP calibration scan -- NO HARDWARE REQUIRED.
%   1. Load the saved scan file
%   2. Dip analysis (simple per-pixel spectral minimum + map)
%   3. FPISpectraAnalysis-based calibration figure (real class, unmodified)
%
% Run this any time after FP_scan_acquire.m has produced a .mat file.
%
% Requires on path: FPISpectraAnalysis.m, Signal Processing Toolbox
% (findpeaks, used inside FPISpectraAnalysis.GetPeaks).

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

fprintf('\n--- Dip analysis summary ---\n');
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
title(sprintf('Center pixel (kx=%d, ky=%d) reflectivity spectrum', kxC, kyC));
legend('Spectrum', 'Dip'); grid on;

subplot(1, 2, 2);
imagesc(scan.x_um, scan.y_um, dipWl_map);
set(gca, 'YDir', 'normal'); axis image; colorbar;
xlabel('x [um]'); ylabel('y [um]'); title('Dip wavelength [nm]');

% Append the dip maps to the loaded file so they don't need recomputing.
try
    save(FILE_TO_ANALYZE, 'dipVal_map', 'dipWl_map', '-append');
catch ME
    warning('Could not append dip maps to %s: %s', FILE_TO_ANALYZE, ME.message);
end

%% ============================== 3. FPISpectraAnalysis (real class) ==============================
% Uses the actual, UNMODIFIED FPISpectraAnalysis.m class. Produces 4 of
% its 6 calibration panels: spectra with min/max peaks and optimal
% wavelength, spectra maxima map, optimal wavelength map, and derivative
% at optimal wavelength map.

smoothSpectraWin    = 10;
smoothDerivativeWin = 20;
smoothCurvatureWin  = 20;

% The class's constructor only avoids a lambdaDim bug when a *trailing*
% size-1 dimension disappears during permute (a quirk of MATLAB's ndims,
% not something the class checks for directly). The original 4D
% acquisition pipeline (nx, ny, avg, nLambda) hit that quirk by accident
% whenever avg==1. Reproduce the same shape here instead of touching the
% class: add back a dummy avg=1 dimension before nLambda, and use
% lambdadim=4 (not 3, which errors on this 3D data).
spectra4D = reshape(spectra, [scan.nx, scan.ny, 1, numel(lambda)]);

fpi = FPISpectraAnalysis('lambda', lambda, 'spectra', spectra4D, 'lambdadim', 4);

fpi.Normalize('spectra');
fpi.Smooth(smoothSpectraWin, 'spectra');
fpi.ComputeDerivative();
fpi.Smooth(smoothDerivativeWin, 'derivative');
fpi.ComputeCurvature();
fpi.Smooth(smoothCurvatureWin, 'curvature');

fpi.GetPeaks();   % MinPeakProminence=0.02, MinPeakHeight=0.15 are hardcoded
                   % inside GetPeaks() itself -- edit them directly in
                   % FPISpectraAnalysis.m to tune peak-detection sensitivity.

fpi.ComputeOptimalWavelengthFromPeaks('maxderivative');

stats = fpi.GetFPICharacteristics();
fpi.ShowFPICharacteristics();

% fpi.ComputeAOTFFrequencies();   % only meaningful if the AOTF calibration
%                                   % file at the hardcoded path in the class
%                                   % constructor actually loaded
% fMHz = [mean(fpi.aotfFrequencies(:)) * 1e-6, std(fpi.aotfFrequencies(:)) * 1e-6];
% fprintf('\n--- FPI analysis summary ---\n');
% fprintf('Mean AOTF frequency: %.3f MHz, std: %.3f MHz\n', fMHz(1), fMHz(2));
% fprintf('Mean optimal wavelength: %.3f nm, std: %.3f nm\n', ...
%     mean(fpi.optimalWavelength(:)) * 1e9, std(fpi.optimalWavelength(:)) * 1e9);

% --- Custom 4-panel figure ---
figure('Name', 'FPI calibration results', 'Position', [50 50 1000 800]);

npx = size(fpi.spectra);
npx = prod(npx(2:end));

subplot(2, 2, 1); hold on;
imagesc(1:npx, fpi.lambda*1e9, fpi.spectra(:,:));
xlabel('Pixel #'); ylabel('Wavelength [nm]');
colormap(gca, 'gray'); colorbar;

for i = 1:npx
    n_min_peaks = length(fpi.peaks.min(i).locs);
    if n_min_peaks
        hm = plot3(i*ones(n_min_peaks,1), fpi.lambda(fpi.peaks.min(i).locs)*1e9, ...
            10*ones(n_min_peaks,1), '.r', 'MarkerSize', 5);
    end
    n_max_peaks = length(fpi.peaks.max(i).locs);
    if n_max_peaks
        hp = plot3(i*ones(n_max_peaks,1), fpi.lambda(fpi.peaks.max(i).locs)*1e9, ...
            10*ones(n_max_peaks,1), 'g.', 'MarkerSize', 5);
    end
    ho = plot3(i, fpi.optimalWavelength(i)*1e9, 10, 'm.', 'MarkerSize', 5);
end
try
    legend([hm hp ho], {'Detected min peaks', 'Detected max peaks', 'Extracted optimal wavelength'});
catch
end
xlim([1 npx]); ylim([min(fpi.lambda) max(fpi.lambda)]*1e9);
set(gca, 'CLim', [0 1]);
title('Normalized spectra'); grid on; box on;

subplot(2, 2, 2);
imagesc(squeeze(fpi.spectraMax));
colorbar; xlabel('x [px]'); ylabel('y [px]');
title('Spectra maxima');

subplot(2, 2, 3);
imagesc(fpi.optimalWavelength*1e9);
colormap(gca, 'turbo'); colorbar;
xlabel('x [px]'); ylabel('y [px]');
title('Extracted optimal wavelength [nm]');

subplot(2, 2, 4);
imagesc(fpi.optimalDerivative);
colormap(gca, 'hot'); colorbar;
xlabel('x [px]'); ylabel('y [px]');
title({'Normalized spectra derivative', 'at optimal wavelength [a.u.]'});

try
    saveas(gcf, fullfile(folder, 'calib_analysis_FPISpectraAnalysis.png'));
catch
end

fprintf('\nAnalysis complete.\n');

%%
%%
fpi = FPISpectraAnalysis('lambda',lambda,'spectra',spectra,'lambdadim',5);

fpi.Normalize('spectra');
fpi.Smooth(20,'spectra');
fpi.ComputeDerivative();
fpi.Smooth(20,'derivative');

fpi.GetPeaks();
fpi.ComputeContrast();
stats = fpi.GetFPICharacteristics();

fpi.ComputeOptimalWavelengthFromPeaks('minderivative');

%%
hf_spec = fpi.ShowSpectra(1,'type','imagesc','showPeaks',true);
hf_spec.CurrentAxes.Title.String = 'Normalized spectra vs time';
hf_spec.CurrentAxes.XLabel.String = 'Time [s]';

t = datetime(times);
t = t - t(1);


%%
peaks.wvl = nan(size(t));
peaks.wid = nan(size(t));

for i=1:length(t)
    if ~isempty(fpi.peaks.min(i).locs)
        if length(fpi.peaks.min(i).locs) > 1
        else
            peaks.wvl(i) = lambda(fpi.peaks.min(i).locs);
            peaks.wid(i) = fpi.peaks.min(i).width;
        end
    end
end

figure(2),clf
    plot(t,peaks.wvl*1e9)
    grid on, box on
    xlabel('Time'),ylabel('Wavelength [nm]')
    title('Peak wavelength vs time')