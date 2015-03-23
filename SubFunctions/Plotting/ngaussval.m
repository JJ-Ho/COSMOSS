function y = ngaussval(x,fwhm)
% FUNCTION y = ngaussval(x,fwhm)
%
% returns the value of a  gaussian with full width at half max fwhm
%   - normalized such that max intensity is 1
%   - integrated peak area is 

a = log(2)/(fwhm/2)^2;

y = exp(-a*x.^2);

end
