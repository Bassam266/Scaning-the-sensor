%% FP_animate_spectra.m
% Animate the reflectivity spectrum at every measurement point of a saved
% FP calibration scan, and save the animation as a GIF.
%
%   1. Load the saved scan file
%   2. Loop over each measurement point (kx, ky) and plot its spectrum
%   3. Save each plotted frame into a GIF
%
% Requires on path: none beyond base MATLAB.

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

lambda_nm = lambda(:) * 1e9;

fprintf('Loaded: %s\n', FILE_TO_ANALYZE);
fprintf('spectra size: %s, lambda size: %s\n', ...
    mat2str(size(spectra)), mat2str(size(lambda)));

%% ============================== 2. SETTINGS ==============================

GIF_FILE   = fullfile(fpath, 'spectra_animation.gif');
FRAME_TIME = 0.1;   % seconds per frame

yMin = min(spectra(:), [], 'omitnan');
yMax = max(spectra(:), [], 'omitnan');

%% ============================== 3. PLOT EACH MEASUREMENT POINT & BUILD GIF ==============================

fig = figure('Name', 'Spectra animation', 'Position', [100 100 700 450]);

isFirstFrame = true;
pt = 0;

for ky = 1:scan.ny
    for kx = 1:scan.nx
        pt = pt + 1;

        pixSpec = squeeze(spectra(kx, ky, :));

        plot(lambda_nm, pixSpec, 'LineWidth', 1.5);
        xlabel('Wavelength [nm]');
        ylabel('Reflectivity');
        ylim([yMin, yMax]);
        grid on;
        title(sprintf('Point %d/%d  (kx=%d, ky=%d)', pt, scan.Npts, kx, ky));

        drawnow;

        % --- capture frame and append to GIF ---
        frame = getframe(fig);
        [imind, cm] = rgb2ind(frame2im(frame), 256);

        if isFirstFrame
            imwrite(imind, cm, GIF_FILE, 'gif', 'Loopcount', inf, 'DelayTime', FRAME_TIME);
            isFirstFrame = false;
        else
            imwrite(imind, cm, GIF_FILE, 'gif', 'WriteMode', 'append', 'DelayTime', FRAME_TIME);
        end
    end
end

close(fig);

fprintf('\nGIF saved to: %s\n', GIF_FILE);
