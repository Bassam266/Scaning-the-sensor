%% FP_calib_peaks_analyze.m
% Analyze a saved FP calibration scan using the FPISpectraAnalysis class:
% normalize/smooth the spectra, extract peaks, compute FPI
% characteristics, and plot the detected minimum-peak wavelength and
% width at every measurement point.
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

%% ============================== 2. FPI ANALYSIS ==============================
% spectra is (nx, ny, nLambda), so wavelength is dimension 3.

fpi = FPISpectraAnalysis('lambda', lambda, 'spectra', spectra, 'lambdadim', 3);

fpi.Normalize('spectra');
fpi.Smooth(20, 'spectra');
fpi.ComputeDerivative();
fpi.Smooth(20, 'derivative');

fpi.GetPeaks();
fpi.ComputeContrast();
stats = fpi.GetFPICharacteristics();

fpi.ComputeOptimalWavelengthFromPeaks('minderivative');

%% ============================== 3. SPECTRA + PEAKS OVERVIEW ==============================

hf_spec = fpi.ShowSpectra(1, 'type', 'imagesc', 'showPeaks', true);
hf_spec.CurrentAxes.Title.String = 'Normalized spectra vs measurement point';
hf_spec.CurrentAxes.XLabel.String = 'Measurement point #';

% No real timestamps are saved for a spatial calibration scan, so use
% the flattened measurement-point index (1:Npts) in place of time.
% obj.peaks.min is linearly indexed the same way as the imagesc pixel
% axis above: index i maps to (kx, ky) in column-major order over
% (scan.nx, scan.ny).
pt = (1:scan.Npts)';

%% ============================== 4. PEAK WAVELENGTH & WIDTH PER POINT ==============================

peaks.wvl = nan(size(pt));
peaks.wid = nan(size(pt));

for i = 1:length(pt)
    if ~isempty(fpi.peaks.min(i).locs)
        if length(fpi.peaks.min(i).locs) > 1
            % Ambiguous: more than one minimum peak found at this point, skip it.
        else
            peaks.wvl(i) = lambda(fpi.peaks.min(i).locs);
            peaks.wid(i) = fpi.peaks.min(i).width;
        end
    end
end

figure(2), clf
plot(pt, peaks.wvl * 1e9)
grid on, box on
xlabel('Measurement point #'), ylabel('Wavelength [nm]')
title('Peak wavelength vs measurement point')
