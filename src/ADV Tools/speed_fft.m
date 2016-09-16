%Brian Polagye
%July 13, 2013

%Description: code to run spectral analysis of velocity data (modified from
%acoustics code)

%Inputs:
%   - y: time domain series to transform to fourier space
%   - t: time stamps associated with y (seconds)
%   - win_pts: number of points in each window (power of two for best
%   performance)
%   - overlap: fractional overlap between windows (0.5 recommended)

%Outputs:
%   - Pyy: power spectral density (units of variance)
%   - f: frequencies associated with Pyy
%   - checksum: confirm variance preserving correction applied correctly

%Change Log:
%  - 4/17/11, BLP - Buffer can create zero entries in first and last column,
%                   prune both
%  - 4/17/13, BLP - Modifications for data processing - primarily
%                   streamlining

function [Pyy, f, checksum] = speed_fft(y, t, win_pts, overlap)

%buffer input time series
y = buffer(y,win_pts,win_pts*overlap);
t = buffer(t,win_pts,win_pts*overlap);

y = y(:,2:end-1);
t = t(:,2:end-1);

%calculate frequency bands
dt = t(2,1)-t(1,1);                     %sampling interval (s)
fs = 1/dt;                              %sampling frequency (Hz)

nfft = size(y,1);                       %length of transform - windows are defined as power of 2
f = fs/2*linspace(0,1,nfft/2)';         %frequency bands - 0:1:Nyquist frequency

%detrend time series with a linear mean
y_detrend = detrend(y,'linear');

%apply hamming window
%note: a rectangular window does a better job of resolving frequency components, but there is
%   high contamination of the spectrum at non-tidal frequencies
%   (Data Analysis Methods in Physical Oceanography p. 448)
wts = repmat(hamming(size(y,1),'periodic'),1,size(y,2));
y_wts = y_detrend.*wts;

S = fft(y_wts,nfft);
%first value is sum of all values (mean) - f = 0 - should be small
%first 1:nfft/2 are positive frequencies up to Nyquist
%second nfft/2:nfft are negative frequencies (complex conjugates of positive frequencies)

%single-sided power spectra (Data Analysis Methods in Physical Oceanography (p. 421), w/ modification
% 2 * complex magnitude(S) squared/(number of points*sampling frequency)
Pyy = 2 * abs(S(1:nfft/2,:)).^2/(nfft*fs);
Pyy(1,:) = Pyy(1,:)/2;      %adjust power in 0 frequency band (mean - should be close to zero)
Pyy(end,:) = Pyy(end,:)/2;  %adjust power in Nyquist band

%adjust power spectral density to satisfy Parseval's theorem E = var(s)
rms_y = (mean(y_detrend.^2)).^0.5;

E = trapz(f,Pyy);      %total energy under spectra                 
corr = E./rms_y.^2;    %empirical correcetion factor

Pyy=Pyy./repmat(corr,size(Pyy,1),1);    %corrected spectral densities

checksum = sum(trapz(f,Pyy)-rms_y.^2);

Pyy = mean(Pyy,2);  %merge all windows to calculate spectrum

end