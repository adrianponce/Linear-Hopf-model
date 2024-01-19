
function [FC,Cov,Cov_t,lags,PowSpect,freq]=StochSim_HopfNet(C,a,g,wo,sigma,T,nTrials,dt)
%
% stochastic network of N hopf nodes
% each node evolves as follows:
%
% dz/dt = (a+iwo)z - z|z|^2 + network_interactions + noise 
%
% network_interactions: g*Cjk(zk-zj)
% 
% Calculates the FC and the power spectral density using stochastic simulations
%
%
% Inputs:
%   - C  : connectivity matrix (N-by-N)
%   - a  : bifurcation parameters for each node (N-dim. vector)
%   - g  : global coupling (scalar)
%   - wo : intrinsic angular frequencies for each node (N-dim. vector)
%   - sigma : noise amplitude (scalar)
%   - varargin (optiona) : lags of the lagged-covariance (default : 0:10 seconds)
%   - T : simulation time after transient dynamics (in seconds)
%   - nTrials : nb. of repetitions to calculate the network statistics
%   - dt : integration step (default: 0.001 seconds)
%
% Outputs:
%   - FC : functional connectivity of real(z)
%   - Cov : covariance matrix of real(z)
%   - Cov_t : lagged covariance of real(z) (N x N x Nlags)
%   - lags : lags of the lagged covariance
%   - PowSpect : power spectral density (PSD) of real(z)
%   - freq : frequencies of the PSD
%
% Adrián Ponce-Alvarez 06-07-2022
%--------------------------------------------------------------------------


% number of nodes:
N = size(C,1);

if nargin<8
dt = 0.001; %in sec
end

tspan = 0:dt:T;
L = length(tspan);

% downsampling:
ds = 20;
Tds=length(0:ds*dt:T);

% number of lags:
nn = 500;

% prepare paramaters for calculating the spectral density
Fs = 1/(dt*ds); % sampling freq

freq = 0:Fs/Tds:Fs/2;

nfreqs=length(freq);
PowSpect=zeros(nfreqs,N);

% Simulation -------------------------------------------------------------
Cov = zeros(N);
Cov_t = zeros(N,N,2*nn+1);

tic
    
    for tr=1:nTrials
        
    fprintf('sim. trial: %g over %g \n',tr,nTrials)   
    Z = zeros(N,Tds);
    z = 0.0001*rand(N,1)+1i*0.0001*rand(N,1);

                %Transient dynamics:
                for t=1:30000
                    r = repmat(z,1,N);
                    z_diffs = r.'-r;
                    u = g*sum(C.*z_diffs,2);
                    z = z + dt*( z.*(a+1i*wo - abs(z).^2) ) ... 
                        + dt*u ...
                        + sqrt(dt)*sigma*randn(N,1) + 1i*sqrt(dt)*sigma*randn(N,1);              
                end
                
                % Stationary regime:
                c=1;
                for t=2:L
                    r = repmat(z,1,N);
                    z_diffs = r.'-r;
                    u = g*sum(C.*z_diffs,2);
                    z = z + dt*( z.*(a+1i*wo - abs(z).^2) ) ... 
                        + dt*u ...
                        + sqrt(dt)*sigma*randn(N,1) + 1i*sqrt(dt)*sigma*randn(N,1);              
                    
                    % down-sampling:
                    if mod(t,ds)==0 
                       Z(:,c)=z; 
                       c=c+1;
                    end
                end       
            
    % covariance of the x-variables (x=real(z)):            
    X=real(Z)';
    Cov = Cov + cov(X)/nTrials;
    
    % power spectrum of real(z) 
     for n=1:N
      xdft = fft(X(:,n));
      ii = floor(Tds/2+1);
      xdft = xdft(1:ii);
      psdx = 2*(1/(Fs*Tds)) * abs(xdft).^2;
%       [psdx,freq] = pspectrum(X(:,n),Fs);
      PowSpect(:,n) = PowSpect(:,n) + psdx/nTrials;
     end
     
     
     % Lagged covariance:
     CovLagged = zeros(N,N,2*nn+1);
     for i=1:N
         for j=1:N
         [clag, lags] = xcov(X(:,i),X(:,j),nn);
         CovLagged(i,j,:) = clag;
         end
     end
     
     Cov_t = Cov_t + CovLagged/nTrials;
    
    end

lags = lags*ds*dt;
Cov_t = Cov_t/Tds;  % normalized lagged-covariance resulting from xcov.m

% correlation matrix:
FC = corrcov(Cov);
comp_time = toc/60;    
fprintf('finished after: %g min \n',comp_time)    
