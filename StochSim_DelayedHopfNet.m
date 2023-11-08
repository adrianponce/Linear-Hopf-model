
function [FC,Cov,Cov_t,lags,PowSpect,freq]=StochSim_DelayedHopfNet(C,D,v,a,g,wo,sigma,T,nTrials,dt)
%
% stochastic network of N hopf nodes
% each node evolves as follows:
%
% dz/dt = (a+iwo)z - z|z|^2 + network_interactions + noise 
%
% network_interactions: g*Cjk*[ zk(t-dkj)-zj(t) ] (dkj: transmission delay)
% 
% Calculates the FC and the power spectral density using stochastic simulations
%
% Inputs:
%  - C : connectivity matrix (N-by-N)
%  - D : delays matrix (N-by-N) in seconds
%  - a : bifurcation parameters for each node (N-dim. vector)
%  - g : global coupling (scalar)
%  - wo : intrinsic angular frequencies for each node (N-dim. vector)
%  - sigma : noise amplitude (scalar)
%  - T : simulation time after transient dynamics (in seconds)
%  - nTrials : nb. of repetitions to calculate the network statistics
%  - dt : integration step (default: 0.001 seconds)
%
% Outputs:
%  - FC : functional connectivity of real(z)
%  - Cov : covariance matrix of real(z)
%  - Cov_t : lagged covariance of real(z) (N x N x Nlags)
%  - lags : lags of the lagged covariance
%  - PowSpect : power spectral density (PSD) of real(z)
%  - freq : frequencies of the PSD
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

% number of nodes:
N = size(C,1);
% node strength
B = sum(C,2);   


if nargin<10
dt = 0.001; %in sec
end

% delays
Delays = round(D/(v*dt));   % Delays matrix in timesteps.
maxD=max(max(Delays)); % get the maximal delay


% add time to reach stationary regime:
Ttran = 500; % in sec

tspan = 0:dt:(T+Ttran);
L = length(tspan);

Ltran = length(0:dt:Ttran);

% downsampling:
ds = 10;
Tds=length(0:ds*dt:(T+Ttran));
Tdstran = length( 0:dt*ds:Ttran );

% number of lags:
nn = 500;

% prepare paramaters for calculating the spectral density
Fs = 1/(dt); % sampling freq
TT = (L-Ltran);
freq = 0:Fs/TT:Fs/2;
%[~,freq] = pspectrum(zeros(1,TT),Fs);
% nfreqs=length(freq);
% PowSpect=zeros(nfreqs,N);


% % prepare paramaters for calculating the spectral density
% Fs = 1/(dt*ds); % sampling freq
% TT = (Tds-Tdstran);
% freq = 0:Fs/TT:Fs/2;
nfreqs=length(freq);
PowSpect=zeros(nfreqs,N);

% Simulation -------------------------------------------------------------
Cov = zeros(N);
Cov_t = zeros(N,N,2*nn+1);


tic
    
    for tr=1:nTrials
        
    fprintf('sim. trial: %g over %g \n',tr,nTrials)   
    Z = zeros(N,Tds);
    z = 0.0001*rand(N,L)+1i*0.0001*rand(N,L);
    
                c=1;
                for t=maxD+1:L-1 
                    
                    u=zeros(N,1);
                    for i=1:N
                        ui=0;
                        for j=1:N
                            ij=Delays(i,j); 
                            ui = ui + g*C(i,j)*z(j,t-ij);                           
                        end
                        u(i)=ui;
                    end
                    
                    z(:,t+1) = z(:,t) + dt*( z(:,t).*(a+1i*wo - abs(z(:,t)).^2 -g*B) ) ... 
                        + dt*u ...
                        + sqrt(dt)*sigma*randn(N,1) + 1i*sqrt(dt)*sigma*randn(N,1);              

                    if mod(t,ds)==0
                       Z(:,c)=z(:,t); 
                       c=c+1;
                    end
                end       

    % remove transient dynamics:
    z(:,1:Ltran) = [];
    Z(:,1:Tdstran) = [];
    
                
    % covariance of the x-variables (x=real(z)): 
    x = real(z)';
    X = real(Z)';
    Cov = Cov + cov(x)/nTrials;
    
    % power spectrum of real(z) 
     for n=1:N
      xdft = fft(x(:,n));
      ii = floor(TT/2+1);
      xdft = xdft(1:ii);
      psdx = 2*(1/(Fs*TT)) * abs(xdft).^2;
%       [psdx,freq] = pspectrum(x(:,n),Fs);
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
Cov_t = Cov_t/(Tds-Tdstran);  % normalized lagged-covariance

% correlation matrix:
FC = corrcov(Cov);
comp_time = toc/60;    
fprintf('finished after: %g min \n',comp_time)    
