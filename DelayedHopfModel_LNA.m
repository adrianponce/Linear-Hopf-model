function [FC,Cov,PSD,freqs,Csp,A]=DelayedHopfModel_LNA(C,D,v,a,g,wo,sigma,varargin)
%
% stochastic network of N hopf nodes
% each node evolves as follows:
%
% dz/dt = (a+iwo)z - z|z|^2 + noise + interactions from the rest of the network
%
% network_interactions: g*Cjk*[ zk(t-dkj)-zj(t) ] (dkj: transmission delay)
% 
% Calculates the FC, the lagged-covariance,the power spectral density, 
% and the cross-spectrum using the linear noise approximation
% (less accurate as the system approaches the bifurcation)
%
% Inputs:
%   - C  : connectivity matrix (N-by-N)
%   - a  : bifurcation parameters for each node (N-dim. vector)
%   - g  : global coupling (scalar)
%   - wo : intrinsic angular frequencies for each node (N-dim. vector)
%   - sigma : noise amplitude (scalar)
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

Del = [D/v zeros(N);zeros(N) D/v]; % in sec


 % Jacobian matrix:
    
    s = sum(C,2);   % node strength
    B = diag( s );
    
    
    
    % no-delays part:
    Axx = diag(a)*eye(N) - g*B;
    Ayy = Axx;
    Axy = -diag(wo);
    Ayx = diag(wo);
    Anod = [Axx Axy; Ayx Ayy];
    
    % delays part:
    Axx_d = g*C;
    Ayy_d = Axx_d;
    Axy_d = zeros(N);
    Ayx_d = zeros(N);
    A_d = [Axx_d Axy_d; Ayx_d Ayy_d];
    
    A = Anod + A_d;
    
    
    % input noise covariance:
    %Qn = sigma^2*eye(2*N);
    

% Power spectrum:
% wo needs to be in Hz!
%--------------------------------------------------------------------------
if nargin <8
df=0.05;
freqs=0:df:50;
else
freqs = varargin{1}; 
df = freqs(2)-freqs(1);
end
numF=length(freqs);

S   = zeros(numF,N);
Csp = zeros(N,N,numF);
I = eye(2*N);

     for k=1:numF
            f   = freqs(k);
            Jaco = Anod + A_d.*exp(+1i*f*(2*pi)*(Del));
            J    = Jaco + 1i*f*(2*pi)*eye(2*N);
            Jt   = J';
            X    = J\I/Jt; %cross-spectrum
            Y   = X(1:N,1:N);
            Csp(:,:,k)=Y;  %cross-spectrum    
            S(k,:) = diag(Y); % spectral density
            % this would be equivalent:
%             Q = inv(J);
%             for i=1:N
%                 S(k,i)=sum(abs(Q(i,:)).^2);
%             end     
     end
            
% psd:
PSD = 2*sigma^2*S;
PSD = real(PSD); % remove null imaginary terms 0.000i

% cross-spectrum:
Csp = 2*sigma^2*Csp;
    
% Covariance = integral of the real part of CS:
Cov = sum(df*real(Csp),3);
FC = corrcov(Cov);
    
 
 return
