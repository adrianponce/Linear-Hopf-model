% This code compares the predictions of the linear approximation against 
% the statistics obtained using stochastic simulations of the nonlinear delayed-coupled model. 
% The coupling matrix was given by the human dMRI connectome from HCP, with N nodes. 
% The interaction delays between nodes were approximated using the Euclidean distance between brain regions 
% divided by a transmission velocity v. We used the distances from the HCP data.
% The model parameters a and ω were drawn from normal distributions N(a0,Δa) and N(ω0,Δω), 
% respectively, with means a0 and ω0, and standard deviations Δa and Δω. 
% We simulated the system for T seconds after letting it reach the stationary regime 
% and we used nTrials realizations of the system with different random initial conditions. 
%
% A Ponce-Alvarez 18-01-2024
%

% load the anatomical connectivity of the brain:
%---------------------------------------------------

load Connectome20.mat % C: structural conn.; D: distances between ROIs [mm]
N = size(C,1);
D = D/1000; % in meters
v = .07; % transmission velocity in m/s 


% Model parameters:
%-----------------------------------------------------

mu_a = -1;
s_a  = 0.3;
a = mu_a + s_a*randn(N,1); % bifurcation parameters
g = 3; % Global coupling

sigma = 0.0002; % noise amplitude
mu_w = 1; % mean frequency
s_w  = .2; % frequency SD
wo = (mu_w + s_w*randn(N,1))*(2*pi);
if any(wo<=0)
    disp('wo < 0 !')
end

% Stochastic simulation:
%-------------------------------------------------------------------------
T  = 600;  % in seconds
nTrials = 100; % number of repetitions of the stochastic simulations
dt = 0.005; % integration step

tic
[FC,Cov,PSD,freq]=StochSim_DelayedHopfNet(C,D,v,a,g,wo,sigma,T,nTrials,dt); 
%   - FC  : correlation matrix of real(z)
%   - Cov : covariance matrix of real(z)
%   - PSD : power spectral density of real(z)
%   - freq : frequencies of the PSD
comp_time = toc/60;    


% Analytical -------------------------------------------------------------
tic
[FC_lna,Cov_lna,PSD_lna,freqs_lna]=DelayedHopfModel_LNA(C,D,v,a,g,wo,sigma,freq);
%   - FC_lna  : correlation matrix of real(z)
%   - Cov_lna : covariance matrix of real(z)
%   - PSD_lna : power spectral density of real(z)
%   - freq_lna : frequencies of the PSD_lna
comp_time_lna = toc/60;    


% Figures:
%--------------------------------------------------------------------------

figure

plot(Cov(:),Cov_lna(:),'k.','markersize',9)
grid on
hold on
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
box on

