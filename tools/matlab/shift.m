function z = shift(x,n,varargin)
%shift - shift x by n samples. Fractions allowed

% Author: Samuel Prager
% University of Southern California
% email: sprager@usc.edu
% Created: 2017/07/29 00:03:16; Last Revised: 2017/07/29 00:03:16

%------------- BEGIN CODE --------------
mode = 'real'; % mode='complex' for complex output even if input is purely real
if (nargin>2)
    mode = varargin{1};
end
if (round(n)==n)
    if (size(x,1)>size(x,2))
        z = circshift(x,[n,0]);
    else
        z = circshift(x,[0,n]);
    end
else
    N = numel(x);
%     t = linspace(-1/2,1/2,N);
    t = [-1/2:1/N:1/2-1/N];
    if (size(x,1)>size(x,2))
        t = t(:);
    end
    maxx = max(abs(x));
    xfft = fftshift(fft(x));
    xfft = xfft.*exp(-1i*2*pi*n*t);
    z = ifft(ifftshift(xfft));

    if(isreal(x) && strcmp(mode,'real'))
        z = real(z);
        scale = sqrt((sum(abs(x).^2))/sum(z.^2));
        z = scale*z;
    end
%     z = (maxx/max(abs(z)))*z;
end

end
%------------- END OF CODE --------------