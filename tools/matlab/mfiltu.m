function [d,l] = mfiltu(x,y,upfac,varargin)
% FUNCTION_NAME - matched filter and upsample. Uses fft/ifft
% Syntax:  [d,l] = mfiltu(x,y,upfac,varargin)
%
% Inputs:
%    x - sample vector
%    y - matched filter samples (taps)
%    upfac - upsample factor (numel(d) = upfac*max(numel(x),numel(y))
%    scalemode - 'none' (default): no output scaling
%                'same': output scaled by upfac to match input scale
%
% Outputs:
%    d - matched filter output. Output length = input length * upfac
%    l - sample lag vector (centered at zero)
%
% Example: 
%    [d,l] = mfiltu(x,y,1)
%    [d,l] = mfiltu(x,y,2,'none')
%    [d,l] = mfiltu(x,y,2,'same')
%
% Other m-files required: fft, ifft
% Subfunctions: none
% MAT-files required: none
%
% See also: xcorr

% Author: Samuel Prager
% University of Southern California
% email: sprager@usc.edu
% Created: 2017/11/07 17:31:18; Last Revised: 2017/11/07 17:31:18

%------------- BEGIN CODE --------------
scalemode = 'none';
if(nargin>3)
    scalemode = varargin{1};
end
% length = max(numel(x),numel(y));
% nzero = (upfac-1)*length/2;
% dfft = fftshift(fft(x,length).*conj(fft(y,length)));
% dfft = [zeros(1,nzero),dfft,zeros(1,nzero)];
% d = fftshift(ifft(ifftshift(dfft)));
% l = [-numel(d)/2:1:(numel(d)/2-1)];
length = max(numel(x),numel(y));
nzero1 = round((upfac-1)*ceil(length/2));
nzero2 = round((upfac-1)*floor(length/2));
dfft = fftshift(fft(x,length).*conj(fft(y,length)));
dfft = [zeros(1,nzero1),dfft,zeros(1,nzero2)];
d = fftshift(ifft(ifftshift(dfft)));
scale = upfac^2;
if(strcmp(scalemode,'same'))
    scale = upfac;
end
d = scale*d;
l = [-numel(d)/2:1:(numel(d)/2-1)];
%------------- END OF CODE --------------
