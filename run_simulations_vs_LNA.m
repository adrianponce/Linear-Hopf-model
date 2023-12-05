% This code compares the predictions of the linear approximation against 
% the statistics obtained using stochastic simulations of the nonlinear model. 
% The coupling matrix was given by the human dMRI connectome from HCP, with N=250 nodes. 
% The model parameters a and ω were drawn from normal distributions N(a0,Δa) and N(ω0,Δω), 
% respectively, with means a0 and ω0, and standard deviations Δa and Δω. 
% We simulated the system for T seconds after letting it reach the stationary regime 
% and we used nTrials realizations of the system with different random initial conditions. 

% load the anatomical connectivity of the brain:
%---------------------------------------------------

load Connectome250.mat
N = size(C,1);

% Model parameters:
%-----------------------------------------------------
mu_a = -1;
s_a  = 0.3;
a = mu_a + s_a*randn(N,1); % bifurcation parameters
g = 3; % Global coupling

sigma = 0.01; % noise amplitude
wo = (1+.2*randn(N,1))*(2*pi); % angular frequencies


% Stochastic simulation:
%-------------------------------------------------------------------------
T  = 600;  % in seconds
nTrials = 100; % number of repetitions of the stochastic simulations

tic
[FC,Cov,Cov_t,lags,PowSpect,freq]=StochSim_HopfNet(C,a,g,wo,sigma,T,nTrials);
%   - FC  : correlation matrix of real(z)
%   - Cov : covariance matrix of real(z)
%   - Ct  : lagged covariance of real(z)
%   - PSD : power spectral density of real(z)
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
% PSD of six example ROIs:
plot(freq,PowSpect(:,1:6))
hold on
col = lines(6);
col = max(col-.2,0);
for i=1:6
plot(freqs_lna,PSD_lna(:,i),'linewidth',2,'color',col(i,:))
end
set(gca,'xlim',[0 3],'xtick',0:.5:3,'linewidth',1)
xlabel('frequency \nu/\nu_{0} (Hz)','fontsize',12)
ylabel('PSD','fontsize',14)


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

xlabel('\tau')
ylabel('C(\tau)')
text(.7,.8,sprintf('node #%g',rs(n)),'fontsize',9,'units','normalized')
end


figure
% Cross-correlation of the most connected pair of ROIs
[~,ind] = max(C(:));
[i,j]=ind2sub(N,ind);

acf_sim = squeeze(Cov_t(i,j,:));
acf_LNA = squeeze(Ct_lna(i,j,:));
plot(lags,acf_sim,'k-','linewidth',2)
hold on
plot(lags_lna,acf_LNA,'r:','linewidth',2)
set(gca,'xlim',[0 3],'fontsize',9)

xlabel('\tau')
ylabel('C_{\itjk\rm}(\tau)')
text(.7,.8,sprintf('nodes #%g and #%g',i,j),'fontsize',9,'units','normalized')


