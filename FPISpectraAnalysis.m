classdef FPISpectraAnalysis < handle
    % Class for scanFP calibration data analysis acquired using the
    % Yokogawa spectrometer
    
    properties (GetAccess = public, SetAccess = private)
        lambda % Optional wavelength vector (must match the corresponding spectra size)
        dlambda = 1; % Optional lambda unit (for multiplying peak widths and everything)
        spectra % N dims spectra matrix with wavelength as the first dimensions
        derivative % First order derivative of the spectra matrix
        curvature % Second order derivative of the spectra matrix
        optimalWavelength % Optimal wavelength matrix
        optimalWavelengthIdx % Optimal wavelength matrix
        optimalWavelengthMap % Optimal wavelength binary map (true if peak has been found, false otherwise)
        optimalDerivative % Values of the derivative at the extracted optimal wavelengths
        wavelenghtAtOptimalPeak % Wavelength at the minima of peaks selected as optimal wavelengths
        normCoeff % Coefficient to normalize the PA signals with (probably obsolete)
        spectraMax % Maximum values of the rawSpectra matrix
        spectraEnergy % integration of each spectra from the rawSpectra matrix
        stats

        aotfFrequencies % Frequencies to drive the AOTF
        aotfCalib_lambda2frequency = @(lambda) ((1973.5547 - lambda*1e9)/11.7783*1e6); % AOTF calibration function
        ase_fit = @(lambda) (0.956*exp(-((lambda-1027e-9)/6.442e-9).^2)+0.0007678);
        
        % Peaks storage sturcture
        peaks = struct('max',struct('val',{},'locs',{},'width',{},'prom',{}),...
                       'min',struct('val',{},'locs',{},'width',{},'prom',{}));
                   
        contrast % Cell array containing contrast of each peaks over the n-dimensions of the spectra
    end
    
    properties (GetAccess = public, SetAccess = protected)
        rawSpectra % Raw measured spectra
    end
    
    properties (GetAccess = private, SetAccess = private)
        lambdaDim = 1; % Spectra wavelength dimension, default is assumed to be 1
        spectraDims % Dimensions of the spectra
                
        isSpectraNormalized = false;
        isSpectraSmoothed = false;
        isSpectraDerivated = false;
        isSpectraCurvated = false;
        
        isDerivativeNormalized = false;
        isDerivativeSmoothed = false;
        
        isCurvatureNormalized = false;
        isCurvatureSmoothed = false;
        
        isOptimalWavelengthComputed = false;
        isGetPeaksComputed = false;
        isContrastComputed = false;
        isNormCoeffComputed = false;
    end
    
    methods
        function obj = FPISpectraAnalysis(varargin)
            % Load the spectra matrix and the optional wavelength vector into the object
            % 
            % Optional input arguments couple:
            %   * 'lambda': should be followed by a wavelength vector
            %   * 'spectra': should be followed by a spectra matrix
            %   * 'lambdaDim': should be followed by an integer indicating
            %   the dimension number of the wavelengths in the spectra
            %   matrix (default is 1)
                
            global mysys

            if ~isempty(varargin) %nargin>1
                for karg = 1:2:length(varargin) %karg=2:2:nargin
                    switch lower(varargin{karg})
                        case 'lambda'
                            obj.lambda = varargin{karg+1};
                        case 'lambdadim'
                            obj.lambdaDim = varargin{karg+1};
                        case 'spectra'
                            obj.spectra = varargin{karg+1};
                            obj.spectraDims = ndims(obj.spectra);
                        otherwise
                            error(['Unknown option ' varargin{karg} ' !'])
                    end
                end
            end
            
            if isempty(obj.spectra)
                warning('Empty spectra property!')
            else
                % Permute the spectra matrix if needed
                if obj.lambdaDim~=1
                    if (obj.lambdaDim <= obj.spectraDims)
                        otherDims = 1:obj.spectraDims;
                        otherDims(otherDims==obj.lambdaDim) = [];
                        obj.spectra = permute(obj.spectra,[obj.lambdaDim, otherDims]);

                        % Adjust obj.spectraDims if last dimension length after
                        % permute is 1 and thus disappears
                        permutedDims = ndims(obj.spectra);
                        if permutedDims~=obj.spectraDims
                            obj.spectraDims = permutedDims;
                        else
                            error('Invalid input argument lambdaDim!')
                        end
                    end
                end
                
                % Check the obj.lambda property validity
                obj.CheckLambda();
                obj.rawSpectra = obj.spectra;
                obj.spectraMax = max(obj.rawSpectra,[],1);
                obj.spectraEnergy = sum(obj.rawSpectra,1);
                
                % Load the aotf calibration file
%                 cal = load('D:\Coding\ChLabGit\pa_rig\aotf\aotfCalib')
% %                 cal = load([mysys.paths.local_rig filesep 'aotf' filesep 'aotfCalib' ])
%                 obj.aotfCalib_lambda2frequency = cal.aotfCalib_lambda2frequency;
            end
            
        end
        
        function CheckLambda(obj)
            % Check if obj.lambda dimension is valid and if it is not empty
            
            if ~isempty(obj.lambda)
                if ~isvector(obj.lambda)
                    error('Please enter a valid lambda vector!')
                else
                    if length(obj.lambda) ~= size(obj.spectra,1)
                        error('Please enter a valid lambda vector!')
                    else
                        obj.dlambda = abs(diff(obj.lambda(1:2)));
                    end
                end
            end
        end
        
        function Normalize(obj,varargin)
            % Normalize the spectra-related matrices by its lambda maximum value
            % 
            % Optional input arguments:
            %   * 'spectra': normalize the spectra matrix
            %   * 'derivative': normalize the derivative matrix
            %   * 'curvature': normalize the curvature matrix
            %
            % If no optional argument is passed, the spectra matrix is
            % normalized by default.
            
            if isempty(obj.spectra)
                error('Empty spectra property!')
            end
            
            if ~isempty(varargin)
                for karg = 1:length(varargin)
                    switch lower(varargin{karg})
                        case 'spectra' % Normalize the spectra
                            maxVals = max(obj.spectra,[],1);
                            obj.spectra = obj.spectra./maxVals;
%                             obj.spectra = obj.spectra ./obj.ase_fit(obj.lambda);
                            obj.isSpectraNormalized = 1;
                        case 'derivative' % Normalize the derivative
                            maxVals = max(obj.derivative,[],1);
                            obj.derivative = obj.derivative./maxVals;
                            obj.isDerivativeNormalized = 1;
                        case 'curvature' % Normalize the curvature
                            maxVals = max(obj.curvature,[],1);
                            obj.derivative = obj.curvature./maxVals;
                            obj.isCurvatureNormalized = 1;
                        otherwise
                            error(['Unknown option ' varargin{karg} ' !'])
                    end
                end
            else % By default, normalize the spectra
                maxVals = max(obj.spectra,[],1);
                obj.spectra = obj.spectra./maxVals;
                obj.isSpectraNormalized = 1;
            end

        end
        
        function Smooth(obj,nMovMean,varargin)
            % Smoothes the spectra-related matrices over the wavelength
            % dimension using moving average function (Matlab's movmean())
            % 
            % Madatory input argument:
            %   * nMovMean: number of points of the moving average
            %
            % Optional input arguments:
            %   * 'spectra': smooth the spectra matrix
            %   * 'derivative': smooth the derivative matrix
            %   * 'curvature': smooth the curvature matrix
            %
            % If no optional argument is passed, the spectra matrix is
            % smoothed by default.
            
            if isempty(obj.spectra)
                error('Empty spectra property!')
            end
            
            if ~isempty(varargin)
                for karg = 1:length(varargin)
                    switch lower(varargin{karg})
                        case 'spectra' % Normalize the spectra
                            obj.spectra = movmean(obj.spectra,nMovMean);
                            obj.isSpectraSmoothed = 1;
                        case 'derivative' % Normalize the derivative
                            obj.derivative = movmean(obj.derivative,nMovMean);
                            obj.isDerivativeSmoothed = 1;
                        case 'curvature' % Normalize the curvature
                            obj.curvature = movmean(obj.curvature,nMovMean);
                            obj.isCurvatureSmoothed = 1;
                        otherwise
                            error(['Unknown option ' varargin{karg} ' !'])
                    end
                end
            else % By default, normalize the spectra
                obj.spectra = movmean(obj.spectra,nMovMean);
                obj.isSpectraSmoothed = 1;
            end
            
        end
        
        function ComputeDerivative(obj)
            % Computes the derivative of the spectra matrix, with respect
            % to the wavelengths dimension.
            %
            % For size compatibility reasons, a nan vector is added to
            % derivative matrix to keep the same dimensions as the spectra
            % matrix.
            
            if isempty(obj.spectra)
                error('Empty spectra property!')
            end
            
            % Compute derivative
            obj.derivative = diff(obj.spectra,1,1);
            
            % Create missing vector
            missingDims = size(obj.spectra);
            missingDims(1) = 1;
            missing = nan(missingDims);
            
            % Add missing vector to match initial spectra matrix dimensions
            obj.derivative = cat(1,missing,obj.derivative);
            obj.isSpectraDerivated = 1;
        end
        
        function ComputeCurvature(obj)
            % Computes the curvature of the spectra matrix, that is the derivative 
            % of the derivative matrix, with respect to the wavelengths dimension.
            %
            % For size compatibility reasons, a nan vector is added to
            % curvature matrix to keep the same dimensions as the spectra
            % matrix.
            
            if isempty(obj.derivative) || obj.isSpectraDerivated~=1
                error('Empty derivative property, please compute it first!')
            end
            
            % Compute derivative
            obj.curvature = diff(obj.derivative,1,1);
            
            % Create missing vector
            missingDims = size(obj.derivative);
            missingDims(1) = 1;
            missing = nan(missingDims);
            
            % Add missing vector to match initial spectra matrix dimensions
            obj.curvature = cat(1,missing,obj.curvature);
            obj.isSpectraCurvated = 1;
        end
        
        
        function ComputeOptimalWavelength(obj,optimalType,opts)
            % Computes the optimal wavelength to bias the interferometer.
            % This computation is based on the derivative maximum values
            % only.
            % To be adjusted with the linear range and the spectral power
            % if we find a proper figure of merit.
            % 
            % Madatory input argument:
            %   * optimalType: string that changes how the optimal
            %   wavelength is computed. Possibles values are the following:
            %       - 'maxabsderivative': returns the wavelengths corresponding to max(abs(obj.derivative))
            %       - 'maxderivative': returns the wavelengths corresponding to max(obj.derivative)
            %       - 'minderivative': returns the wavelengths corresponding to min(obj.derivative)
            %       - 'midlinearrange': requires threshold linThresh
            %       - 'midminmax'
            % Optional input argument:
            %   * lambdaRAnge: wavelength range
            %   * linThresh: threshold for linear range
            
            arguments
                obj
                optimalType (1,:) char {mustBeMember(optimalType,{'maxabsderivative',...
                    'maxderivative','minderivative','midminmax','midlinearrange'})}
                opts.lambdaRange (1,2) {mustBeNumeric}
                opts.linThresh (1,1) {mustBeNumeric} = 0.8
            end
            
            if isempty(obj.derivative)
                error('Please compute the derivative of your spectra!')
            end
            
            lambda = obj.lambda;
            spectra = obj.spectra; % Store spectra matrix
            derivative = obj.derivative; % Store derivative matrix
            
            if isfield(opts,'lambdaRange')
                [idx,~,~] = find(obj.lambda>=opts.lambdaRange(1) & obj.lambda<=opts.lambdaRange(2));
                lambda = lambda(idx);
                
                if ~isempty(idx)
                    sz = size(derivative);
                    derivative = derivative(idx,:); % Crop wavelengths and auto reshape
                    derivative = reshape(derivative,[length(idx) sz(2:end)]); % Reshape to new wavelength dimensions with other initial dimensions
                    spectra = spectra(idx,:);
                    spectra = reshape(spectra,[length(idx) sz(2:end)]);
                end
            end
            
            switch lower(optimalType)
                case 'maxabsderivative' % Returns the wavelengths corresponding to max(abs(obj.derivative))
                    [obj.optimalDerivative,slopeIndex] = max(abs(derivative),[],1,'omitnan');
                case 'maxderivative' % Returns the wavelengths corresponding to max(obj.derivative)
                    [obj.optimalDerivative,slopeIndex] = max(derivative,[],1,'omitnan');
                case 'minderivative' % Returns the wavelengths corresponding to min(obj.derivative)
                    [obj.optimalDerivative,slopeIndex] = min(derivative,[],1,'omitnan');
                
                case 'midlinearrange' % to be vectorized
                    derivative = derivative(:,:);
                    for k = 1:size(derivative,2)
                        [idx,~,~] = find(derivative(:,k)>=opts.linThresh*max(derivative(:,k)));
                        idxStart(k) = idx(1);
                        idxEnd(k) = idx(end);
                    end
                    slopeIndex = idxStart + floor(0.5*(idxEnd-idxStart));
                    slopeIndex = reshape(slopeIndex,sqrt(size(derivative,2)),sqrt(size(derivative,2)));
                                    
                case 'midminmax'
                    derivative = derivative(:,:);
                    for k = 1:size(derivative,2)
                        
                         % cut derivative from after min
                        [~,slopeIndexMax] = max(derivative(:,k),[],1,'omitnan');
                        [~,slopeIndexMin] = min(derivative(:,k),[],1,'omitnan');
                        
                        % find index when derivative is crossing 0 right after min and right after max
                        [idx,~,~] = find(derivative(slopeIndexMin:end,k)>=0,1,'first');
                        idx1(k) = idx + slopeIndexMin-1;
                                                
                        [idx,~,~] = find(derivative(slopeIndexMax:end,k)<0,1,'first');
                        idx2(k) = idx + slopeIndexMax-1;
                    end
                    
                    slopeIndex = idx1 + floor(0.5*(idx2-idx1));
                    slopeIndex = reshape(slopeIndex,sqrt(size(derivative,2)),sqrt(size(derivative,2)));

                otherwise
                    error(['Unknown optimal wavelength type ' optimalType '!']);
            end
            
            slopeIndex = squeeze(slopeIndex);
            
            obj.optimalWavelength = lambda(slopeIndex);
            obj.optimalWavelengthIdx = slopeIndex;
            obj.isOptimalWavelengthComputed = 1;
            obj.optimalDerivative = obj.derivative(obj.optimalWavelengthIdx);
            
            [~,peakIndex] = min(abs(spectra),[],1,'omitnan');
            obj.wavelenghtAtOptimalPeak = lambda(peakIndex);
        end
        
        
        function ComputeOptimalWavelengthFromPeaks(obj,optimalType,varargin)
            % Computes the optimal wavelength to bias the interferometer.
            % This computation is based on the derivative maximum values
            % around each peaks. The getPeaks method should be prealably
            % invoked.
            % To be adjusted with the linear range and the spectral power
            % if we find a proper figure of merit.
            % 
            % Madatory input argument:
            %   * optimatType: string that changes how the optimal
            %   wavelength is computed. Possibles values are the following:
            %       - 'maxabsderivative': returns the wavelengths corresponding to max(abs(obj.derivative))
            %       - 'maxderivative': returns the wavelengths corresponding to max(obj.derivative)
            %       - 'minderivative': returns the wavelengths corresponding to min(obj.derivative)
            %
            % Optional input argument:
            %   * varargin{1}: wavelength range
            
            if isempty(obj.derivative)
                error('Please compute the derivative of your spectra!')
            end
            
            if ~obj.isGetPeaksComputed
                error('Please invoke the GetPeaks() method first')
            end
            
            
            %% Start by reshaping the derivative by the user input wavelengths if need be
            lambda = obj.lambda;
            spectra = obj.spectra; % Store spectra matrix
            derivative = obj.derivative; % Store derivative matrix
            
            narginchk(1,3);
            if ~isempty(varargin)
                lambdaRange = varargin{1};
                [idx,~,~] = find(obj.lambda>=lambdaRange(1) & obj.lambda<=lambdaRange(2));
                lambda = lambda(idx);
                
                if ~isempty(idx)
                    sz = size(derivative);
                    derivative = derivative(idx,:); % Crop wavelengths and auto reshape
                    derivative = reshape(derivative,[length(idx) sz(2:end)]); % Reshape to new wavelength dimensions with other initial dimensions
                    spectra = spectra(idx,:);
                    spectra = reshape(spectra,[length(idx) sz(2:end)]);
                end
            end
            
            
            %% Restrict the optimal derivative range depending on the peak width
            npx = size(obj.spectra);
            npx = prod(npx(2:end)); % Number of pixels to analyze
            lambda_idx = cell(npx,1); % Indices of the wavelengths around each peak
            
            for i=1:npx
                n_min_peaks = length(obj.peaks.min(i).locs);
                if n_min_peaks
                    if n_min_peaks >= 2 % Get min peak index
                        [~,idx_pk] = min(obj.peaks.min(i).val);
                    else
                        idx_pk = 1;
                    end
                    
                    try
                        [lambda_idx{i},~,~] = find( abs(obj.lambda - obj.lambda(obj.peaks.min(i).locs(idx_pk))) < 3*obj.peaks.min(i).width(idx_pk) );
                    catch
                        disp(['Skipped pixel number ' num2str(i)])
                    end
                end
            end
                
            
            %% Compute the optimal derivative based on user input
            sz_spec = size(obj.spectra);
            slopeIndex = nan(sz_spec(2:end));
            obj.optimalWavelength = nan(sz_spec(2:end));
            obj.optimalWavelengthIdx = nan(sz_spec(2:end));
            obj.optimalWavelengthMap = nan(sz_spec(2:end));
            obj.optimalDerivative = nan(sz_spec(2:end));
            obj.wavelenghtAtOptimalPeak = nan(sz_spec(2:end));
            
            switch lower(optimalType)
                case 'maxabsderivative' % Returns the wavelengths corresponding to max(abs(obj.derivative))
                    for i=1:npx
                        if ~isempty(lambda_idx{i})
                            [obj.optimalDerivative(i),sidx] = max(abs(derivative(lambda_idx{i},i)),[],1,'omitnan');
                            slopeIndex(i) = lambda_idx{i}(sidx);
                        end
                    end
                    
                case 'maxderivative' % Returns the wavelengths corresponding to max(obj.derivative)
                    for i=1:npx
                        if ~isempty(lambda_idx{i})
                            [obj.optimalDerivative(i),sidx] = max(derivative(lambda_idx{i},i),[],1,'omitnan');
                            slopeIndex(i) = lambda_idx{i}(sidx);
                        end
                    end
                    
                case 'minderivative' % Returns the wavelengths corresponding to min(obj.derivative)
                    for i=1:npx
                        if ~isempty(lambda_idx{i})
                            [obj.optimalDerivative(i),sidx] = min(derivative(lambda_idx{i},i),[],1,'omitnan');
                            slopeIndex(i) = lambda_idx{i}(sidx);
                        end
                    end
                    
                otherwise
                    error(['Unknown optimal wavelength type ' optimalType '!']);
            end
            
            valid_idx = ~isnan(slopeIndex);
            invalid_idx = isnan(slopeIndex);
            obj.optimalWavelength(valid_idx) = lambda(slopeIndex(valid_idx));
            obj.optimalWavelength(invalid_idx) = 1028e-9; % Default value (needed by the AOTF)
            obj.optimalWavelengthIdx = slopeIndex;
            obj.optimalWavelengthMap(valid_idx) = true;
            obj.optimalWavelengthMap(invalid_idx) = false;
            obj.isOptimalWavelengthComputed = 1;
            
            for i=1:npx
                if ~isempty(lambda_idx{i})
                    [~,pidx] = min(spectra(lambda_idx{i},i),[],1,'omitnan');
                    obj.wavelenghtAtOptimalPeak(i) = lambda(lambda_idx{i}(pidx));
                end
            end
        end
        
        
        function ComputeAOTFFrequencies(obj,varargin)
            % Computes the frequencies to send to the AOTF from the
            % computed optimal wavelengths.
            % 
            % Optional input arguments:
            %   * varargin couple: 'function' followed by a custom
            %   function handle to convert the wavelengths into
            %   frequencies from calibration data.
            
            if ~obj.isOptimalWavelengthComputed
                error('Please compute the optimal wavelength matrix first!')
            end
            
            functionFlag = 0;
            
            if ~isempty(varargin)
                for karg = 1:2:length(varargin)
                    switch lower(varargin{karg})
                        case 'function' % Use the custom function for conversion
                            functionFlag = 1;
                            func = varargin{karg+1};
                        otherwise
                            error(['Unknown option ' varargin{karg} ' !'])
                    end
                end
            end
            
            if functionFlag
                obj.aotfFrequencies = func(obj.optimalWavelength);
            else
                obj.aotfFrequencies = obj.aotfCalib_lambda2frequency(obj.optimalWavelength);
            end
            % Proper calibration to be added, pass optional anonymous function
        end
        
       
        function GetPeaks(obj)
            % Gets the peaks with respect to the wavelength dimension for the spectra
            % matrix
            
            if ~obj.isSpectraNormalized
                obj.Normalize('spectra');
            end
            
            % Get the initial size of the spectra, lambda should be the first dimension
            sizeSpectra = size(obj.spectra);
            
            % Initalize cell to be all the remaining dimensions to automatically cover the N-dimensional case
            sizeCell = prod(sizeSpectra(2:end));
            initCell = cell(sizeCell,1);
            
            % Initialize a structure with cell fields to contain all peaks informations
            peaksMax = struct('val',initCell,'locs',initCell,'width',initCell,'prom',initCell);
            peaksMin = struct('val',initCell,'locs',initCell,'width',initCell,'prom',initCell); 
            
            % Loop through all dimensions to find all maxima and minima peaks
            for k = 1:sizeCell
                [pks,locs,width,prom] = findpeaks(obj.spectra(:,k),'MinPeakProminence',0.02,'MinPeakHeight',0.15);
                peaksMax(k).val = pks;
                peaksMax(k).locs = locs;
                peaksMax(k).width = width*obj.dlambda;
                peaksMax(k).prom = prom;

                [pks,locs,width,prom] = findpeaks(1-obj.spectra(:,k),'MinPeakProminence',0.02,'MinPeakHeight',0.15);
                peaksMin(k).val = 1-pks;
                peaksMin(k).locs = locs;
                peaksMin(k).width = width*obj.dlambda;
                peaksMin(k).prom = prom;
            end
            
            % Reshape to the original spectra dimensions
            if length(size(obj.spectra)) > 2
                peaksMax = reshape(peaksMax,sizeSpectra(2:end));
                peaksMin = reshape(peaksMin,sizeSpectra(2:end));
            end
            
            % Store resuts into the dedicated object peaks property
            obj.peaks.max = peaksMax;
            obj.peaks.min = peaksMin;
            obj.isGetPeaksComputed = 1;
        end
        
        function ComputeContrast(obj)
            % Computes the contrast of the peaks for any spectral
            % dimension (to be refined).
            %
            % Super difficult to auto;ate perfectly, see the comments
            % below.
            
            if ~obj.isGetPeaksComputed
                error('Please get peaks data first using the obj.GetPeaks() method!')
            else
                if isempty(obj.peaks.max) || isempty(obj.peaks.min)
                    error('No peaks have been found, aborting!')
                end
            end
            
            % Get the initial size of the spectra, lambda should be the first dimension
            sizeSpectra = size(obj.spectra);
            
            % Initalize cell to be all the remaining dimensions to automatically cover the N-dimensional case
            sizeCell = prod(sizeSpectra(2:end));
            C = cell(sizeCell,1);
            
            % Loop through all dimensions to find compute all contrast from peaks
            for k = 1:sizeCell
                % For now it is assumed that a detected maximum peak should
                % necessarily be followed by a minimum peak, and that they
                % should correspond to the same peak from which the
                % contrast is computed. Of course this might not be 
                % generally the case...
                % 
                % Really hard to automate, to be refined, but better than
                % nothing although not perfect
                % 
                % Array crossing and peaks locations could certainly be
                % used it it really works like shit
                
                % Get the current number of min and max peaks values and
                % the array size difference
                minVals = obj.peaks.min(k).val;
                maxVals = obj.peaks.max(k).val;
                
                nmin = length(minVals);
                nmax = length(maxVals);
                dn = nmax-nmin;
                
                if (nmin==0) || (nmax==0)
                    C{k} = nan;
                    continue
                end
                
                % Fill the shortest array with nans (dn==0 means same lengths)
                if dn>0
                    minVals = cat(1,minVals,nan(abs(dn),1));
                elseif dn<0
                    maxVals = cat(1,maxVals,nan(abs(dn),1));
                end
                    
                % Compute the contrast value
                C{k} = (maxVals - minVals)./(maxVals + minVals);
            end
            
            % Reshape to the original spectra dimensions
            if length(size(obj.spectra)) > 2
                C = reshape(C,sizeSpectra(2:end));
            end
            
            % Store resuts into the dedicated object contrast property
            obj.contrast = C;
            obj.isContrastComputed = 1;
        end
        
        function stats = GetFPICharacteristics(obj)
            % Computes the FPI peak characterstics such as: contrast, 
            % width, slope

            % Get statistics about spectra max amplitudes
            stats.spectraAmplitudes.mean = mean(obj.spectraMax, 'all', 'omitnan');
            stats.spectraAmplitudes.std = std(obj.spectraMax, 0, 'all', 'omitnan');
            stats.spectraAmplitudes.max = max(obj.spectraMax, [], 'all', 'omitnan');
            stats.spectraAmplitudes.min = min(obj.spectraMax, [], 'all', 'omitnan');
            
            % Get statistics about optimal wavelengths, validity maps, and derivatives
            stats.optimalWavelengths = [];
            stats.optimalWavelengthsMap = [];
            stats.optimalDerivatives = [];
            stats.wavelengthsAtOptimalPeaks = [];
            if (obj.isSpectraDerivated && obj.isOptimalWavelengthComputed)
                stats.wavelengthsAtOptimalPeaks.mean = mean(obj.wavelenghtAtOptimalPeak, 'all', 'omitnan');
                stats.wavelengthsAtOptimalPeaks.std = std(obj.wavelenghtAtOptimalPeak, 0, 'all', 'omitnan');
                stats.wavelengthsAtOptimalPeaks.max = max(obj.wavelenghtAtOptimalPeak, [], 'all', 'omitnan');
                stats.wavelengthsAtOptimalPeaks.min = min(obj.wavelenghtAtOptimalPeak, [], 'all', 'omitnan');
                
                stats.optimalWavelengths.mean = mean(obj.optimalWavelength, 'all', 'omitnan');
                stats.optimalWavelengths.std = std(obj.optimalWavelength, 0, 'all', 'omitnan');
                stats.optimalWavelengths.max = max(obj.optimalWavelength, [], 'all', 'omitnan');
                stats.optimalWavelengths.min = min(obj.optimalWavelength, [], 'all', 'omitnan');
                
                stats.optimalDerivatives.mean = mean(obj.optimalDerivative, 'all', 'omitnan');
                stats.optimalDerivatives.std = std(obj.optimalDerivative, 0, 'all', 'omitnan');
                stats.optimalDerivatives.max = max(obj.optimalDerivative, [], 'all', 'omitnan');
                stats.optimalDerivatives.min = min(obj.optimalDerivative, [], 'all', 'omitnan');
                
                stats.optimalWavelengthsMap.valid = numel(obj.optimalWavelengthMap == true);
                stats.optimalWavelengthsMap.invalid = numel(obj.optimalWavelengthMap == false);
                stats.optimalWavelengthsMap.total = numel(obj.optimalWavelengthMap);
                stats.optimalWavelengthsMap.ratioValid = stats.optimalWavelengthsMap.valid/stats.optimalWavelengthsMap.total;
                stats.optimalWavelengthsMap.ratioInvalid = stats.optimalWavelengthsMap.invalid/stats.optimalWavelengthsMap.total;
            end

            
            % Get statistics about contrast
            stats.contrast = [];
            if (obj.isContrastComputed && ~isempty(obj.contrast))
                Ctmp = [];
                for k=1:numel(obj.contrast)
                    Ctmp = [Ctmp; obj.contrast{k}];
                end
                
                stats.contrast.mean = mean(Ctmp, 'omitnan');
                stats.contrast.std = std(Ctmp, 0, 1, 'omitnan');
                stats.contrast.max = max(Ctmp, 'all', 'omitnan');
                stats.contrast.min = min(Ctmp, 'all', 'omitnan');
            end
            
            
            % Get statistics about peak widths
            stats.peaks = [];
            if (obj.isGetPeaksComputed && ~isempty(obj.peaks.min))
                Wtmp = [];
                for k=1:numel(obj.peaks.min)
                    Wtmp = [Wtmp; obj.peaks.min(k).width];
                end
                
                stats.peaks.width.mean = mean(Wtmp, 'omitnan');
                stats.peaks.width.std = std(Wtmp, 0, 1, 'omitnan');
                stats.peaks.width.max = max(Wtmp, 'all', 'omitnan');
                stats.peaks.width.min = min(Wtmp, 'all', 'omitnan');
            end
            
            obj.stats = stats;
            
        end
        
        
        function ShowFPICharacteristics(obj)
            if ~isempty(obj.stats)                
                fprintf('\n\n\t Spectra statistics (computed from %d points):', numel(obj.spectra)/size(obj.spectra,1) )
                
                if ~isempty(obj.stats.spectraAmplitudes)
                    fprintf('\n\t\t > Spectra amplitudes [uW]: ');
                    fprintf('mean = %1.1f, std = %1.1f, max = %1.1f, min = %1.1f', ...
                        obj.stats.spectraAmplitudes.mean*1e6, obj.stats.spectraAmplitudes.std*1e6, ...
                        obj.stats.spectraAmplitudes.max*1e6, obj.stats.spectraAmplitudes.min*1e6);
                end
                
                if ~isempty(obj.stats.optimalWavelengthsMap)
                    fprintf('\n\t\t > Percentage of points with a cavity peak: %1.1f', obj.stats.optimalWavelengthsMap.ratioValid*100);
                end
                
                if ~isempty(obj.stats.optimalWavelengths)
                    fprintf('\n\t\t > Optimal wavelengths [nm]: ');
                    fprintf('mean = %1.1f, std = %1.1f, max = %1.1f, min = %1.1f', ...
                        obj.stats.optimalWavelengths.mean*1e9, obj.stats.optimalWavelengths.std*1e9, ...
                        obj.stats.optimalWavelengths.max*1e9, obj.stats.optimalWavelengths.min*1e9);
                end
                
                if ~isempty(obj.stats.wavelengthsAtOptimalPeaks)
                    fprintf('\n\t\t > Wavelengths at optimal peaks [nm]: ');
                    fprintf('mean = %1.1f, std = %1.1f, max = %1.1f, min = %1.1f', ...
                        obj.stats.wavelengthsAtOptimalPeaks.mean*1e9, obj.stats.wavelengthsAtOptimalPeaks.std*1e9, ...
                        obj.stats.wavelengthsAtOptimalPeaks.max*1e9, obj.stats.wavelengthsAtOptimalPeaks.min*1e9);
                end
                
                if ~isempty(obj.stats.optimalDerivatives)
                    fprintf('\n\t\t > Derivatives at optimal wavelengths [a.u.]: ');
                    fprintf('mean = %1.2e, std = %1.2e, max = %1.2e, min = %1.2e', ...
                        obj.stats.optimalDerivatives.mean, obj.stats.optimalDerivatives.std, ...
                        obj.stats.optimalDerivatives.max, obj.stats.optimalDerivatives.min);
                end
                
                fprintf('\n\n');
            else
                warning('No statistics to show yet, compute them first! Aborting.');
            end
            
        end
        
        function GetNormCoeff(obj)
            if obj.isOptimalWavelengthComputed
                obj.normCoeff = nan(size(obj.optimalWavelength,2),size(obj.optimalWavelength,3));
                
                for kx=1:size(obj.optimalWavelengthIdx,2)
                    for ky=1:size(obj.optimalWavelengthIdx,3)
                        obj.normCoeff(kx,ky) = obj.rawSpectra(obj.optimalWavelengthIdx(kx,ky),kx,ky);
                    end
                end
                obj.normCoeff = obj.normCoeff/max(obj.normCoeff,[],'all');
                
                obj.isNormCoeffComputed = 1;
            else
                error('Please compute the optimal wavelengths first.');
            end
        end
        
        
        function hfig = ShowSpectra(obj,fignum,varargin)
            % Show spectra characteristics
            %
            % Mandatory input:
            %   * fignum: figure number to plot on
            % 
            % Optional inputs:
            %   * 'type': string followed by the plot type as string: either 'plot' or 'imagesc'
            %   * 'showpeaks': string followed by a boolean
            %
            % Output:
            %   * hfig: figure handle
            
            plotTypes = {'plot','imagesc'};
            
            % Initialize with default values
            plotType = 'plot';
            showPeaks = false;
            
            % Check for optional arguments
            while ~isempty(varargin)
                switch lower(varargin{1})
                    case 'type'
                        if any(contains(plotTypes,lower(varargin{2})))
                            plotType = lower(varargin{2});
                        else
                            error(['Unknown plot type ' lower(varargin{2}) '.']);
                        end
                    case 'showpeaks'
                        showPeaks = varargin{2} & obj.isGetPeaksComputed;
                    otherwise
                        warning(['Unknown optional argument ' varargin{1} '. Discarded.'])
                end
                varargin(1:2) = [];
            end
            
            
            % Plot according to the selected plot type
            if ~isempty(obj.spectra)
                hfig = figure(fignum); clf, hold on
                
                if strcmpi(plotType,'plot')
                    plot(obj.lambda*1e9,obj.spectra(:,:));
                    xlabel('Wavelength [nm]')
                    ylabel('Normalized spectra [a.u.]')
                    
                elseif strcmpi(plotType,'imagesc')
                    npx = size(obj.spectra);
                    npx = prod(npx(2:end));
                    imagesc(1:npx,obj.lambda*1e9,obj.spectra(:,:));
                    xlabel('Pixel #')
                    ylabel('Wavelength [nm]')
                    colormap gray
                    colorbar
                    
                    if showPeaks % Display computed peak locations
                        for i=1:npx
                            n_min_peaks = length(obj.peaks.min(i).locs);
                            if n_min_peaks
                                x_vec = i*ones(n_min_peaks,1);
                                y_vec = obj.lambda(obj.peaks.min(i).locs)*1e9;
                                z_vec = 10*ones(n_min_peaks,1);
                                hm = plot3(x_vec,y_vec,z_vec,'.r','Markersize',5);
                            end
                            
                            n_max_peaks = length(obj.peaks.max(i).locs);
                            if n_max_peaks
                                x_vec = i*ones(n_max_peaks,1);
                                y_vec = obj.lambda(obj.peaks.max(i).locs)*1e9;
                                z_vec = 10*ones(n_max_peaks,1);
                                hp = plot3(x_vec,y_vec,z_vec,'g.','Markersize',5);
                            end
                            
                            if obj.isOptimalWavelengthComputed
                                ho = plot3(i,obj.optimalWavelength(i)*1e9,10,'m.','Markersize',5);
                            end
                        end
                        
                        if obj.isOptimalWavelengthComputed
                            legend([hm hp ho],{'Detected min peaks','Detected max peaks','Extracted optimal wavelength'})
                        else
                            legend([hm hp],{'Detected min peaks','Detected max peaks'})
                        end
                    end
                    
                    xlim([1 npx])
                    ylim([min(obj.lambda) max(obj.lambda)]*1e9)
                    caxis([0 1])
                    grid on, box on
                else
                    error('Unknown plot type');
                end
                
                title('Normalized spectra')
                grid on, box on
            else
                error(['spectra property is empty']);
            end
            
            
        end
        
        
        
        function hfig = ShowSpectraEnergy(obj,fignum)
            % Show spectra energy
            %
            % Mandatory input:
            %   * fignum: figure number to plot on
            % 
            % Output:
            %   * hfig: figure handle
            
            % Plot according to the selected optima type
            hfig = figure(fignum); clf
            imagesc(squeeze(obj.spectraEnergy))
            title('Spectra energy')
            colorbar, colormap gray
            xlabel('x [px]'),ylabel('y [px]')
            grid on, box on
        end
        
        
        function hfig = ShowSpectraMax(obj,fignum)
            % Show spectra max
            %
            % Mandatory input:
            %   * fignum: figure number to plot on
            % 
            % Output:
            %   * hfig: figure handle
            
            % Plot according to the selected optima type
            hfig = figure(fignum); clf
            imagesc(squeeze(obj.spectraMax))
            title('Spectra max [W]')
            colorbar, colormap gray
            xlabel('x [px]'),ylabel('y [px]')
            grid on, box on
        end

        
        
        function hfig = ShowDerivative(obj,fignum,varargin)
            % Show derivative characteristics
            %
            % Mandatory input:
            %   * fignum: figure number to plot on
            % 
            % Optional inputs:
            %   * 'type': string followed by the plot type as string: either 'plot' or 'imagesc'
            %
            % Output:
            %   * hfig: figure handle
            
            plotTypes = {'plot','imagesc'};
            
            % Initialize with default values
            plotType = 'plot';
            
            % Check for optional arguments
            while ~isempty(varargin)
                switch lower(varargin{1})
                    case 'type'
                        if any(contains(plotTypes,lower(varargin{2})))
                            plotType = lower(varargin{2});
                        else
                            error(['Unknown plot type ' lower(varargin{2}) '.']);
                        end
                    otherwise
                        warning(['Unknown optional argument ' varargin{1} '. Discarded.'])
                end
                varargin(1:2) = [];
            end
            
            
            % Plot according to the selected plot type
            if ~isempty(obj.derivative)
                hfig = figure(fignum); clf
                
                if strcmpi(plotType,'plot')
                    plot(obj.lambda*1e9,obj.derivative(:,:));
                    xlabel('Wavelength [nm]')
                    ylabel('Normalized derivative [a.u.]')
                elseif strcmpi(plotType,'imagesc')
                    npx = size(obj.spectra);
                    npx = prod(npx(2:end));
                    imagesc(1:npx,obj.lambda*1e9,obj.derivative(:,:));
                    xlabel('Pixel #')
                    ylabel('Wavelength [nm]')
                    colorbar
                else
                    error('Unknown plot type');
                end
                
                title('Spectra')
                grid on, box on
            else
                error(['derivative property is empty']);
            end
            
            
        end
        
        
        
        function hfig = ShowCurvature(obj,fignum,varargin)
            % Show curvature characteristics
            %
            % Mandatory input:
            %   * fignum: figure number to plot on
            % 
            % Optional inputs:
            %   * 'type': string followed by the plot type as string: either 'plot' or 'imagesc'
            %
            % Output:
            %   * hfig: figure handle
            
            plotTypes = {'plot','imagesc'};
            
            % Initialize with default values
            plotType = 'plot';
            
            % Check for optional arguments
            while ~isempty(varargin)
                switch lower(varargin{1})
                    case 'type'
                        if any(contains(plotTypes,lower(varargin{2})))
                            plotType = lower(varargin{2});
                        else
                            error(['Unknown plot type ' lower(varargin{2}) '.']);
                        end
                    otherwise
                        warning(['Unknown optional argument ' varargin{1} '. Discarded.'])
                end
                varargin(1:2) = [];
            end
            
            
            % Plot according to the selected plot type
            if ~isempty(obj.curvature)
                hfig = figure(fignum); clf
                
                if strcmpi(plotType,'plot')
                    plot(obj.lambda*1e9,obj.curvature(:,:));
                    xlabel('Wavelength [nm]')
                    ylabel('Normalized curvature [a.u.]')
                elseif strcmpi(plotType,'imagesc')
                    npx = size(obj.spectra);
                    npx = prod(npx(2:end));
                    imagesc(1:npx,obj.lambda*1e9,obj.curvature(:,:));
                    xlabel('Pixel #')
                    ylabel('Wavelength [nm]')
                    colorbar
                else
                    error('Unknown plot type');
                end
                
                title('Curvature')
                grid on, box on
            else
                error(['curvature property is empty']);
            end
            
            
        end
        
        
        function hfig = ShowOptima(obj,fignum,varargin)
            % Show spectra optima
            %
            % Mandatory input:
            %   * fignum: figure number to plot on
            % 
            % Optional inputs:
            %   * 'type': string followed by the optima type as string: either 'wavelength', 'slope' or 'curvature'
            %
            % Output:
            %   * hfig: figure handle
            
            optTypes = {'wavelength','derivative','curvature'};
            
            % Initialize with default values
            optType = 'wavelength';
            
            % Check for optional arguments
            while ~isempty(varargin)
                switch lower(varargin{1})
                    case 'type'
                        if any(contains(optTypes,lower(varargin{2})))
                            optType = lower(varargin{2});
                        else
                            error(['Unknown optima type ' lower(varargin{2}) '.']);
                        end
                    otherwise
                        warning(['Unknown optional argument ' varargin{1} '. Discarded.'])
                end
                varargin(1:2) = [];
            end
            
            
            
            % Plot according to the selected optima type
            if obj.isOptimalWavelengthComputed
                hfig = figure(fignum); clf
                
                if strcmpi(optType,'wavelength')
                    imagesc(squeeze(obj.optimalWavelength)*1e9)
                    title('Optimal wavelength [nm]')
                    
                elseif strcmpi(optType,'derivative')
                    if ~isempty(obj.derivative)
                        imagesc(squeeze(obj.derivative(obj.optimalWavelengthIdx)))
                        title('Derivative at optimal wavelength')
                    else
                        error('Please compute the derivative first.');
                    end
                    
                elseif strcmpi(optType,'curvature')
                    if ~isempty(obj.curvature)
                        imagesc(squeeze(obj.curvature(obj.optimalWavelengthIdx)))
                        title('Curvature at optimal wavelength')
                    else
                        error('Please compute the curvature first.');
                    end
                else
                    error('Unknown optima type');
                end
                
                colorbar
                xlabel('x [px]'),ylabel('y [px]')
                grid on, box on
            else
                error('Please compute the optimal wavelength first.');
            end

        end


        function hfig = ShowCalibrationResults(obj, opts)
            arguments
                obj FPISpectraAnalysis
                opts.fignum (1,1) double {mustBeInteger, mustBePositive, mustBeNonzero} = 1
                opts.newfig (1,1) logical = 0
            end
            
            if opts.newfig
                hfig = figure(); clf
            else
                hfig = figure(opts.fignum); clf
            end
            
            %                 figure()
                hax(1) = subplot(2,3,1); hold on
                    npx = size(obj.spectra);
                    npx = prod(npx(2:end));
                    imagesc(1:npx,obj.lambda*1e9,obj.spectra(:,:));
                    xlabel('Pixel #')
                    ylabel('Wavelength [nm]')
                    colormap gray
                    colorbar
                    
                    for i=1:npx
                        n_min_peaks = length(obj.peaks.min(i).locs);
                        if n_min_peaks
                            x_vec = i*ones(n_min_peaks,1);
                            y_vec = obj.lambda(obj.peaks.min(i).locs)*1e9;
                            z_vec = 10*ones(n_min_peaks,1);
                            hm = plot3(x_vec,y_vec,z_vec,'.r','Markersize',5);
                        end
                        
                        n_max_peaks = length(obj.peaks.max(i).locs);
                        if n_max_peaks
                            x_vec = i*ones(n_max_peaks,1);
                            y_vec = obj.lambda(obj.peaks.max(i).locs)*1e9;
                            z_vec = 10*ones(n_max_peaks,1);
                            hp = plot3(x_vec,y_vec,z_vec,'g.','Markersize',5);
                        end
                        
                        if obj.isOptimalWavelengthComputed
                            ho = plot3(i,obj.optimalWavelength(i)*1e9,10,'m.','Markersize',5);
                        end
                    end
                    
                    try
                        if obj.isOptimalWavelengthComputed
                            legend([hm hp ho],{'Detected min peaks','Detected max peaks','Extracted optimal wavelength'})
                        else
                            legend([hm hp],{'Detected min peaks','Detected max peaks'})
                        end
                    end
                    xlim([1 npx])
                    ylim([min(obj.lambda) max(obj.lambda)]*1e9)
                    caxis([0 1])
                    grid on, box on
                
                    title('Normalized spectra')
                    grid on, box on


%                 figure()
                hax(2) = subplot(2,3,2);
                    imagesc(squeeze(obj.spectraMax))
                    colorbar
                    xlabel('x [px]'),ylabel('y [px]')
                    title('Spectra maxima [uW]')
% figure()
                hax(3) = subplot(2,3,3);
                    imagesc(obj.optimalWavelength*1e9)
                    colormap(hax(3), 'parula')
% colormap parula
                    colorbar
                    xlabel('x [px]'),ylabel('y [px]')
                    title('Extracted optimal wavelength [nm]')
% figure()
                hax(4) = subplot(2,3,4);
                    imagesc(obj.optimalWavelengthMap)
                    colormap(hax(4), 'gray')
% colormap gray
                    colorbar
                    xlabel('x [px]'),ylabel('y [px]')
                    title('Peak validity binary map')

% figure()
                hax(6) = subplot(2,3,6);
                    imagesc(obj.optimalDerivative)
                    colormap(hax(6), 'hot')
%                     colormap hot
                    colorbar
                    xlabel('x [px]'),ylabel('y [px]')
                    title({'Normalized spectra derivative', 'at optimal wavelength [a.u.]'})

        end
        

    end
    
    
    methods (Static)
        function str = GetStrFromStructVal(s, opts)
            % Returns a string from a struct containing values in each field
            %
            % Input: struct
            % Output: string
            
            arguments
                s struct
                opts.delimiter char = ', '
                opts.joiner char = ' = '
                opts.unitstring char = []
                opts.multiplier double = 1
            end
            
            str = [];
            fields = fieldnames(s);
            
            for f = 1:length(fields)
                current = [ fields{f} opts.joiner num2str(s.(fields{f})*opts.multiplier,'%5.2f') opts.unitstring];
                if f ~= length(fields)
                    current = [current opts.delimiter];
                end
                    
                str = [str, current];
            end
        end
    end
end
