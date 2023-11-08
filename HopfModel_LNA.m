function [FC,Cov,Ct,lags,PSD,freqs,Csp,A]=HopfModel_LNA(C,a,g,wo,sigma,varargin)
%
% stochastic network of N hopf nodes
% each node evolves as follows:
%
% dz/dt = (a+iwo)z - z|z|^2 + noise + interactions from the rest of the network
%
% network interactions: g*Cjk(zk-zj)
% 
% Calculates the FC, the lagged-covariance,the power spectral density, 
% and the cross-spectrum using the linear noise approximation
% (if the origin z = 0 is stable)
% (less accurate as the system approaches the bifurcation)
%
% Inputs:
%   - C  : connectivity matrix (N-by-N)
%   - a  : bifurcation parameters for each node (N-dim. vector)
%   - g  : global coupling (scalar)
%   - wo : intrinsic angular frequencies for each node (N-dim. vector)
%   - sigma : noise amplitude (scalar)
%   - varargin (optiona) : lags of the lagged-covariance (default : 0:10 seconds)
%
% Outputs:
%   - FC  : correlation matrix of real(z)
%   - Cov : covariance matrix of real(z)
%   - Ct  : lagged covariance of real(z)
%   - PSD : power spectral density of real(z)
%   - lags: lags used for Ct
%   - freqs: frequencies used for PSD
%   - Csp : cross-spectrum
%   - A : Jacobian matrix
%
% Adrián Ponce-Alvarez 06-07-2022
%--------------------------------------------------------------------------

% ensure that a and wo are column vectors:
if ~iscolumn(a)
    a = transpose(a);
end
if ~iscolumn(wo)
    wo = transpose(wo);
end

N = size(C,1);

 % Jacobian matrix:
    
    s = sum(C,2);   % node strength
    B = diag( s );
    
    Axx = diag(a)*eye(N) - g*B + g*C;
    Ayy = Axx;
    Axy = -diag(wo);
    Ayx = diag(wo);
    
    A = [Axx Axy; Ayx Ayy];
    
    % input noise covariance:
    Qn = sigma^2*eye(2*N);
    
% Check stability of the origin:
[~,d] = eig(A);
d = diag(d);
Remax = max(real(d));
if Remax >= 0
   disp('Warning: the origin is not stable') 
   FC = []; Cov=[]; Ct=[]; lags=[]; PSD=[]; freqs=[]; Csp=[];
   return
end
    

% Covariance equation:
Cv = lyap(A,Qn); % Solves the Lyapunov equation: A*Cv + Cv*A' + Qn = 0
Corr=corrcov(Cv);

FC = Corr(1:N,1:N); % correlation matrix of real(Z)
Cov = Cv(1:N,1:N); % covariance matrix of real(Z)

% Lagged covariance:
%--------------------------------------------------------------------------
if nargin < 6
L=10;
lags=0:.1:L; % default lags to evaluate the lagged-covariances
else
lags = varargin{1}; % if asked, use these lags    
end
num_lags = length(lags);
Ct=zeros(N,N,num_lags);

 for n=1:length(lags)
     t = lags(n);
     Y=expm(A*t)*Cv;
     Ct(:,:,n) = Y(1:N,1:N);
 end


% Power spectrum:
% wo needs to be in radians per second!
%--------------------------------------------------------------------------
if nargin < 7
dfr=0.005;
freqs=0.05:dfr:6; % default frequencies to evaluate the PSDs and cross-spectrum
else
freqs = varargin{2}; % if asked, use these frequencies   
end
numF=length(freqs);

S   = zeros(numF,N);
Csp = zeros(N,N,numF,'single');

     for k=1:numF
            f   = freqs(k);
            J   = A+1i*f*(2*pi)*eye(2*N);
            Q   = inv(J);
            X   = Q*(Q'); 
            Y   = X(1:N,1:N);
            Csp(:,:,k) = single(Y);  %cross-spectrum    
            S(k,:) = diag(Y); % spectral density
            % this would be equivalent:
            %for i=1:N
            %    S(k,i)=sum(abs(Q(i,:)).^2);
            %end     
     end
            
% PSD:
PSD = 2*sigma^2*S;

% cross-spectrum:
Csp = 2*sigma^2*Csp;
 
 return
