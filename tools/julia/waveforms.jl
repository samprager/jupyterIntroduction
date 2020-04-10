module Waveforms

using FFTW

export chirp, mfiltu, zeropad, zeroext

# function chirp(t,bw)
# Description:
#   Generate LFM chirp from time index t with bandwidth bw
function chirp(t,bw)
    return exp.(1im*pi*.5*(bw/t[end])*(t.^2));
end

# function zeropad(x,nzero,mode::String="back")
# Description:
#    Add nzero zeros to vector x.
# Inputs:
#    x - sample vector
#    nzero - number of zeros to add
#    mode - "back" (default): append zeros to back of vector
#           "front": preappend zeros to front of vector
function zeropad(x,nzero,mode::String="back")
    if (mode=="back")
        xpad = [x; zeros(eltype(x),max(0,Int(nzero)))];
    else
        xpad = [zeros(eltype(x),max(0,Int(nzero)));x];
    end
    return xpad;
end

# function zeroext(x,n,mode::String="back")
# Description:
#    Extend vector x to length n by adding zeros
# Inputs:
#    x - sample vector
#    n - desired length of x
#    mode - "back" (default): append zeros to back of vector
#           "front": preappend zeros to front of vector
function zeroext(x,n,mode::String="back")
    if (mode=="back")
        xpad = [x; zeros(eltype(x),max(0,Int(n)-length(x)))];
    else
        xpad = [zeros(eltype(x),max(0,Int(n)-length(x)));x];
    end
    return xpad;
end


# function mfiltu(x,y,upfac=1,scalemode::String="none")
# Description:
#   Matched filter x with y and upsample by factor upfac
# Inputs:
#    x - sample vector
#    y - matched filter samples (taps)
#    upfac - upsample factor (numel(d) = upfac*max(numel(x),numel(y))
#    scalemode - "none" (default): no output scaling
#                "same": output scaled by upfac to match input scale
# Outputs:
#    d - matched filter output. Output length = input length * upfac
#    l - sample lag vector (centered at zero)
#
# Example:
#    d,l = mfiltu(x,y,1)

function mfiltu(x,y,upfac=1,scalemode::String="none")
    len = max(length(x),length(y));
    nzero1 = Int(round((upfac-1)*ceil(len/2)));
    nzero2 = Int(round((upfac-1)*floor(len/2)));
    # xpad = [x; zeros(eltype(x),max(0,len-length(x)))];
    # ypad = [y; zeros(eltype(y),max(0,len-length(y)))];
    xpad = zeroext(x,len);
    ypad = zeroext(y,len);
    dfft = fftshift(fft(xpad).*conj(fft(ypad)));
    dfft = [zeros(eltype(dfft),nzero1);dfft;zeros(eltype(dfft),nzero2)];
    d = fftshift(ifft(ifftshift(dfft)));
    scale = Float64(upfac^2);
    if(scalemode=="same")
        scale = Float64(upfac);
    end
    d = scale*d;
    l = ((-length(d)/2):1:(length(d)/2-1));
    return d,l;
end

end
