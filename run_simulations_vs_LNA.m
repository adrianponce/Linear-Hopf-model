% This code compares the predictions of the linear approximation against 
% the statistics obtained using stochastic simulations of the nonlinear model. 
% The coupling matrix was given by the human dMRI connectome from HCP, with N nodes. 
% The model parameters a and ω were drawn from normal distributions N(a0,Δa) and N(ω0,Δω), 
% respectively, with means a0 and ω0, and standard deviations Δa and Δω. 
% We simulated the system for T seconds after letting it reach the stationary regime 
% and we used nTrials realizations of the system with different random initial conditions. 
%
% A Ponce-Alvarez 18-01-2024
%

% load the anatomical connectivity of the brain:
%---------------------------------------------------

load Connectome20.mat
N = size(C,1);

% Model parameters:
%-----------------------------------------------------
mu_a = -1;
s_a  = 0.3;
a = mu_a + s_a*randn(N,1); % bifurcation parameters
g = 3; % Global coupling

sigma = 0.01; % noise amplitude

mu_w = 1; % mean frequency
s_w  = .2; % frequency SD
wo = (mu_w + s_w*randn(N,1))*(2*pi); % angular frequencies


% Stochastic simulation:
%-------------------------------------------------------------------------
T  = 600;  % in seconds
nTrials = 100; % number of repetitions of the stochastic simulations

tic
[FC,Cov,Cov_t,lags,PSD,freq]=StochSim_HopfNet(C,a,g,wo,sigma,T,nTrials);
%   - FC : functional connectivity of real(z)
%   - Cov : covariance matrix of real(z)
%   - Cov_t : lagged covariance of real(z) (N x N x Nlags)
%   - lags : lags of the lagged covariance
%   - PSD : power spectral density (PSD) of real(z)
%   - freq : frequencies of the PSD
comp_time = toc/60;    


% Analytical -------------------------------------------------------------
testfreqs = freq; % tested frequencies to evaluate the PSDs
tic
[FC_lna,Cov_lna,Ct_lna,lags_lna,PSD_lna,freqs_lna,CS_lna,A]= ...
    HopfModel_LNA(C,a,g,wo,sigma,lags,testfreqs);
%   - FC_lna  : correlation matrix of real(z)
%   - Cov_lna : covariance matrix of real(z)
%   - Ct_lna  : lagged covariance of real(z)
%   - PSD_lna : power spectral density of real(z)
%   - Csp     : cross-spectrum
%   - A       : Jacobian matrix
comp_time_lna = toc/60;    


% Figures:
%--------------------------------------------------------------------------

figure
% Covariances:
plot(Cov(:),Cov_lna(:),'k.','markersize',9)
hold on
grid on
xlabel('Covariances (simulation)')
ylabel('Linear prediction')
ylim = get(gca,'ylim');
plot(ylim,ylim,'k:')


figure
% PSD of six example ROIs (randomly chosen):
col = lines(6);
hold on
rs = randsample(1:N,6);
for i = 1:6
    x = freq/mu_w;
    y = PSD(:,rs(i))/(sigma^2);
    plot(x,y,'color',min(col(i,:)+.15,ones(1,3)))
    x = freqs_lna/mu_w;
    y = PSD_lna(:,rs(i))/(sigma^2);
    plot(freqs_lna,y,'color',max(col(i,:)-.2,zeros(1,3)),'linewidth',2)
end
set(gca,'xlim',[0 5],'xtick',0:1:5,'fontsize',10,'linewidth',1)
xlabel('frequency \nu/\nu_{0}','fontsize',11)
ylabel('PSD:  \phi(\nu)/\sigma^2','fontsize',11)


figure
% Autocorrelation of three example ROIs:
rs = randsample(1:N,3);
for n = 1:3
axes('position',[.12 .75-(n-1)*.3 .8 .2])    

acf_sim = squeeze(Cov_t(rs(n),rs(n),:));
acf_LNA = squeeze(Ct_lna(rs(n),rs(n),:));
plot(lags,acf_sim,'k-','linewidth',2)
hold on
plot(lags_lna,acf_LNA,'r:','linewidth',2)

set(gca,'xlim',[0 3],'fontsize',9)
xlabel('lag \itt\rm [s]')
ylabel('C(\itt\rm)')
end


figure
% Cross-correlation of a pair of ROIs
rs = randsample(1:N,2);
i = rs(1);
j = rs(2);
acf_sim = squeeze(Cov_t(i,j,:));
acf_LNA = squeeze(Ct_lna(i,j,:));
plot(lags,acf_sim,'k-','linewidth',2)
hold on
plot(lags_lna,acf_LNA,'r:','linewidth',2)
set(gca,'xlim',[0 3],'fontsize',9)

xlabel('lag \itt\rm [s]')
ylabel('C_{\itjk\rm}(\itt\rm)')


