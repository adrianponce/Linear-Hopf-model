% This code compares the predictions of the linear approximation against 
% the statistics obtained using stochastic simulations of the nonlinear delayed-coupled model. 
% The coupling matrix was given by the human dMRI connectome from HCP, with N=250 nodes. 
% The interaction delays between nodes were approximated using the Euclidean distance between brain regions 
% divided by a transmission velocity v. We used the distances from the HCP data.
% The model parameters a and ω were drawn from normal distributions N(a0,Δa) and N(ω0,Δω), 
% respectively, with means a0 and ω0, and standard deviations Δa and Δω. 
% We simulated the system for T seconds after letting it reach the stationary regime 
% and we used nTrials realizations of the system with different random initial conditions. 

% load the anatomical connectivity of the brain:
%---------------------------------------------------

load Connectome250.mat % C: structural conn.; D: distances between ROIs [mm]
N = size(C,1);
D = D/1000; % in meters
v = .07; % transmission velocity in m/s 


% Model parameters:
%-----------------------------------------------------

mu_a = -1;
s_a  = 0.5;
a = mu_a + s_a*randn(N,1); % bifurcation parameters
g = 15; % Global coupling

sigma = 0.0002; % noise amplitude
mu_w = 1;
s_w  = .3;
wo = (mu_w + s_w*randn(N,1))*(2*pi);
if any(wo<=0)
    disp('wo < 0 !')
end


% Stochastic simulation:
%-------------------------------------------------------------------------
T  = 600;  % in seconds
nTrials = 100; % number of repetitions of the stochastic simulations

tic
[FC,Cov,Cov_t,lags,PSD,freq]=StochSim_DelayedHopfNet(C,D,v,a,g,wo,sigma,T,nTrials); 
%   - FC  : correlation matrix of real(z)
%   - Cov : covariance matrix of real(z)
%   - Ct  : lagged covariance of real(z)
%   - PSD : power spectral density of real(z)
comp_time = toc/60;    


% Analytical -------------------------------------------------------------
testfreqs = freq; % tested frequencies to evaluate the PSDs
tic
[FC_lna,Cov_lna,PSD_lna,freqs_lna]=DelayedHopfModel_LNA(C,D,v,a,g,wo,sigma,testfreqs);
%   - FC_lna  : correlation matrix of real(z)
%   - Cov_lna : covariance matrix of real(z)
%   - Ct_lna  : lagged covariance of real(z)
%   - PSD_lna : power spectral density of real(z)
comp_time_lna = toc/60;    


% Figures:
%--------------------------------------------------------------------------

figure

plot(Cov(:),Cov_lna(:),'k.','markersize',9)
hold on
ii = find(eye(size(Cov)));
plot(Cov(ii),Cov_lna(ii),'.','color',[.2 .2 1],'markersize',9)
grid on
xlabel('Covariances (simulation)')
ylabel('Linear prediction')
ylim = get(gca,'ylim');
plot(ylim,ylim,'k:')
text(.07,.9,'variances <(\delta\itx_{j}\rm)^2>','color',[.2 .2 1],'fontsize',9,'units','normalized')
text(.07,.82,'covariances <\delta\itx_{j}\rm\delta\itx_{k}\rm>','fontsize',9,'units','normalized')


figure
% PSD of four example ROIs:
col = lines(5);
hold on
for i = 1:4
    y = PSD(:,i);
    plot(freq,y,'color',min(col(i,:)+.15,ones(1,3)))
    y = PSD_lna(:,i);
    plot(freqs_lna,y,'color',max(col(i,:)-.2,zeros(1,3)),'linewidth',2)
end
set(gca,'xlim',[0 6],'xtick',0:1:6,'linewidth',1)
set(gca,'fontsize',10)
xlabel('frequency \nu/\nu_{0}','fontsize',11)
ylabel('PDF','fontsize',11)
box on


