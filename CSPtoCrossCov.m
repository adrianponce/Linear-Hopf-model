function [Cov_t,lags] = CSPtoCrossCov(Csp,f,lags)
%
% Gets cross-covariance from cross-spectrum
% The cross-covariance is the inverse Fourier transform of the
% cross-spectrum.
%
% Inputs:
% - Csp : cross-spectrum, size NxNxNf, where N is the number of nodes and
% Nf is the number of frequencies
% - f : frequencies
% - lags : lags to calculate the cross-covariance
% 
% Outputs:
% - Cov_t_lna : cross-covariance
% - lags_lna : lags to calculate the cross-covariance (only positive)
%
% Adrián Ponce-Alvarez 06-12-2023
%--------------------------------------------------------------------------

N = size(Csp,1);

% ensure that f is column vector:
if ~iscolumn(f)
    f = transpose(f);
end
df = f(2)-f(1);
tt = lags(lags>=0); % take only positive lags
lags = tt;
L = length(tt);
Cov_t = zeros(N,N,L);

for i=1:N
    for j=1:N
    csp = squeeze(Csp(i,j,:));  
        % Inverse Fourier transform:
        for k = 1:L
        t = tt(k);    
        eiwt = exp( 1i*f*(2*pi)*t );
        Cov_t(i,j,k) = sum(df*csp.*eiwt);
        end
    end
end
Cov_t = real(Cov_t);

return
