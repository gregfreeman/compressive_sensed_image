function [m_theta, samples]=bcs_gf(Phi, v, Ms, IdxParent, IdxChildren, wLevel,hyperpara, MCMCpara, plotflag)
%TSWCSmcmc_NoSc: Tree-structured wavelet CS inversion implemented by MCMC, without scaling coefficients.
%USAGE: [theta, samples]=TSWCSmcmc_NoSc(Phi, v, Ms, IdxParent, IdxChildren, hyperpara, MCMCpara, plotflag) 
%INPUT (number in [] means default value):  
%   Phi: N x M, CS projection matrix
%   v: N x 1, CS observation
%   Ms: 1 x L, number of wavelet coefficients for each decomposition level (scaling coeff. are excluded) 
%       (L represents the number of decomposition levels)   
%   IdxParent: M x 1, parent index. IdxParent(i)=j means the parent of the ith wavelet coefficient is wavelet coefficient j 
%   IdxChildren: (M-Ms(L)) x nChildren, children index. IdxChildren(i,:)=[j(1),...,j(nChildren)] means
%                the children of the ith wavelet coefficient are wavelet coefficients j(1),...,j(nChildren). Since wavelet
%                coefficients at level L (finest level) do not have any children, they are not stored.  
%   hyperpara: MCMC hyperparameters
%       hyperpara.a: scalar, hyperparameter 1 for noise precision [1e-6]
%       hyperpara.b: scalar, hyperparameter 2 for noise precision [1e-6]
%       hyperpara.c: scalar, hyperparameter 1 for non-zero coefficients precision [1e-6] 
%       hyperpara.d: scalar, hyperparameter 2 for non-zero coefficients precision [1e-6] 
%       hyperpara.er: scalar, hyperparameter 1 for weights of root wavelet coefficients [0.9*Ms(1)] 
%       hyperpara.fr: scalar, hyperparameter 2 for weights of root wavelet coefficients [0.1*Ms(1)]  
%       hyperpara.e0: 1 x L, hyperparameters 1 for weights of wavelet coefficients with zero parent for each 
%                     wavelet level. e0(1) is useless, convenient for coding. [ [0, 1/M, ..., 1/M].*Ms ]
%       hyperpara.f0: 1 x L, hyperparameters 2 for weights of wavelet coefficients with zero parent for each 
%                     wavelet level. f0(1) is useless, convenient for coding. [ [0, 1-1/M, ..., 1-1/M].*Ms ]
%       hyperpara.e1: 1 x L, hyperparameters 1 for weights of wavelet coefficients with nonzero parent for each 
%                     wavelet level. e1(1) is useless, convenient for coding. [ [0, 0.5, ..., 0.5].*Ms ]
%       hyperpara.f1: 1 x L, hyperparameters 2 for weights of wavelet coefficients with nonzero parent for each 
%                     wavelet level. f1(1) is useless, convenient for coding. [ [0, 0.5, ..., 0.5].*Ms ]
%   MCMCpara: MCMC parameters
%       MCMCpara.nBurnin: scalar, number of burnin iterations [200]
%       MCMCpara.nCollect: scalar, number of collected samples [100]
%       MCMCpara.thinFactor: scalar, samples are collected every MCMCpara.thinFactor iterations [1]  
%   plotflag: indicator, one means plotting result for each iteration, zero means no plot [0] 
%OUTPUT:
%   theta: M x 1, mean of reconstructed wavelet coefficients (scaling coefficients are excluded)
%   samples: nCollect x 1 cell, collected samples for posterior distribution
%       samples(i).theta: M x 1, theta sample for the ith collection
%       samples(i).trans: 2 x 2 x L, transition (weights) sample for each level for the ith collection.   
%                  samples(i).trans(1,:,1) is null, used for scaling coeff. weights in TSWCSmcmc_NoSc, 
%                  samples(i).trans(2,:,1) is root wavelet coeff. weights for [zero, nonzero], 
%                  samples(i).trans(1,:,k) is weights of [zero, nonzero] for wavelet coeff. at level k (k>1) with zero parent, 
%                  samples(i).trans(2,:,k) is weights of [zero, nonzero] for wavelet coeff. at level k (k>1) with nonzero parent, 

%--------------------------------------------------------------------------
% References:
% L.He and L.Carin, "Exploiting Structure in Wavelet-Based Bayesian Compressive Sensing" (2008)  
%
% Lihan He, ECE, Duke University
% Created: Dec. 19, 2008
% Last change: Mar. 3, 2009, allowing [] for input arguments (using default value)
%--------------------------------------------------------------------------

% ---------------------
% check input arguments
% ---------------------

if nargin<9, plotflag=0; end
if nargin<8
    MCMCpara.nBurnin=200;
    MCMCpara.nCollect=100;
    MCMCpara.thinFactor=1;
end
if nargin<7
    % Number of wavelet coefficients
    M=size(Phi,2);
    % Number of decomposition level
    L=length(Ms);
    % hyperparameters
    hyperpara.a=1e-6;
    hyperpara.b=1e-6;
    hyperpara.c=1e-6;
    hyperpara.d=1e-6;
    hyperpara.er=0.9*Ms(1);
    hyperpara.fr=0.1*Ms(1);
    hyperpara.e0=[0, ones(1, L-1)/M].*Ms;
    hyperpara.f0=[0, 1-ones(1, L-1)/M].*Ms;
    hyperpara.e1=[0, ones(1, L-1)*0.5].*Ms;
    hyperpara.f1=[0, ones(1, L-1)*0.5].*Ms;
end

if isempty(hyperpara)
    % Number of wavelet coefficients
    M=size(Phi,2);
    % Number of decomposition level
    L=length(Ms);
    % hyperparameters
    hyperpara.a=1e-6;
    hyperpara.b=1e-6;
    hyperpara.c=1e-6;
    hyperpara.d=1e-6;
    hyperpara.er=0.9*Ms(1);
    hyperpara.fr=0.1*Ms(1);
    hyperpara.e0=[0, ones(1, L-1)/M].*Ms;
    hyperpara.f0=[0, 1-ones(1, L-1)/M].*Ms;
    hyperpara.e1=[0, ones(1, L-1)*0.5].*Ms;
    hyperpara.f1=[0, ones(1, L-1)*0.5].*Ms;
end
if isempty(MCMCpara)
    MCMCpara.nBurnin=200;
    MCMCpara.nCollect=100;
    MCMCpara.thinFactor=1;
end
if isempty(plotflag), plotflag=0; end

a=hyperpara.a;
b=hyperpara.b;
c=hyperpara.c;
d=hyperpara.d;
er=hyperpara.er;
fr=hyperpara.fr;
e0=hyperpara.e0;
f0=hyperpara.f0;
e1=hyperpara.e1;
f1=hyperpara.f1;

std_v=std(v);
v=v/std_v;

% -------------------
% Data specifications
% -------------------

% Number of CS measurements
N=length(v);       
% Number of wavelet coefficients
M=size(Phi,2);
% Number of decomposizion level
L=length(Ms);

% --------------
% Initialization
% --------------

% theta: M x 1, estimated wavelet coefficients
theta=zeros(M,1);
pi_tilde=zeros(M,1);
% alpha: 1 x L, coefficient precision, coefficients at each wavelet level share one common alpha
alpha=ones(1,L);

% PI: M x 1, mixing weight for each wavelet coefficient
PI=er/(er+fr)*ones(M,1);
for s=2:L
    PI(wLevel==s)=e0(s)/(e0(s)+f0(s))*ones(Ms(s),1);
end

% alpha0: scalar, noise precision
alpha0=1/(std(v)^2/1e2);

%---------------
% precomputation
%---------------

% Phi_{i}^{T}Phi_{i} for all i, M x 1 vector
% PhiTPhi=sum(Phi.*Phi,1)';
PhiTPhi=size(Phi,1)./size(Phi,2); % approx? given random ?

% --------------
% Gibbs Sampling
% --------------

for iter=1:(MCMCpara.nBurnin+MCMCpara.nCollect)

    % (1) theta -- sequentially drawn
    % \tilde{alpha}_{i}, M x 1
    alpha_tilde=alpha(wLevel)'+alpha0*PhiTPhi;
    % i=1
%     v_res=Phi*theta-v;
%     v_tilde=v-A+Phi(:,1)*theta(1);
%     mu_tilde=(alpha0*Phi(:,1)'*v_tilde)/alpha_tilde(1);

    mu=alpha0.*(1./alpha(wLevel)').*(Phi'*v);
    for s=1:L
        ratio=sqrt(alpha(wLevel)'./alpha_tilde).*exp(0.5*alpha_tilde.*mu.*mu).*PI./(1-PI);
        infratio=isinf(ratio);
        pi_tilde(infratio)=1;
        pi_tilde(~infratio)=ratio(~infratio)./(ratio(~infratio)+1);
        u=rand(M,1);
        nonzero=u<pi_tilde;
        theta(~nonzero & wLevel==s)=0;
        theta(nonzero & wLevel==s)=normrnd(mu(nonzero& wLevel==s), sqrt(1./alpha_tilde(nonzero& wLevel==s)));
    end

    % (2) alpha
    for s=1:L
        nonzero_s=nonzero & (wLevel==s);
        sum_idx=sum(nonzero_s);
        if sum_idx>0
            alpha(s)=gamrnd(c+0.5*sum_idx, 1/(d+0.5*sum(theta(nonzero_s).^2,1)));
        end
    end
        
    % (3) PI
%     indi=zeros(M,1); indi(z)=1;     % nonzero coefficient indicator
    Trans=zeros(2,2,L);             % store current transition matrix
    % for root level
    wLevel_s=wLevel==1; %idxLevelStart(1):idxLevelEnd(1);
    Nnon0=sum(wLevel_s&nonzero);
    beta_sample=betarnd(er+Nnon0,fr+sum(wLevel_s)-Nnon0);
    PI(wLevel_s)=beta_sample;
    Trans(2,:,1)=[1-beta_sample, beta_sample];
    wLevel_s_parents=wLevel_s;
    % for other levels
    for s=2:L
        wLevel_s=wLevel==s;
        nLevel_s=sum(wLevel_s);
        wLevel_s_has_nonzero_parent=false(M,1);
        wLevel_s_has_nonzero_parent(wLevel_s)=nonzero(IdxParent(wLevel_s));
        % for parent=0
        Nnon0=sum(wLevel_s_parents & ~nonzero);
        beta_sample=betarnd(e0(s)+Nnon0, f0(s)+nLevel_s-Nnon0);
        PI(nLevel_s & ~wLevel_s_has_nonzero_parent)=beta_sample;
        Trans(1,:,s)=[1-beta_sample, beta_sample];
        % for parent~=0
        Nnon0=sum(wLevel_s_parents & nonzero);
        beta_sample=betarnd(e1(s)+Nnon0, f1(s)+nLevel_s-Nnon0);
        PI(wLevel_s_has_nonzero_parent)=beta_sample;
        Trans(2,:,s)=[1-beta_sample, beta_sample];
        wLevel_s_parents=wLevel_s;
    end
    PI(PI==1)=1-eps;
    PI(PI==0)=eps;
    
    % (4) alpha0
    res=v-Phi(:,nonzero)*theta(nonzero);
    alpha0=gamrnd(a+N/2, 1/(b+res'*res/2));
    
    % Collect samples
    if iter>MCMCpara.nBurnin && (mod(iter-MCMCpara.nBurnin,MCMCpara.thinFactor)==1 || MCMCpara.thinFactor==1)
        i = ceil((iter-MCMCpara.nBurnin)/MCMCpara.thinFactor);
        samples(i).theta = theta*std_v;
        samples(i).trans = Trans;
    end

    fprintf('Iteration %d.\n',iter);
    
    % plot
    if plotflag
        figure(2000), 
        subplot(3,1,1), plot(theta*std_v,'.'); axis([0,M+1,1.1*min(theta*std_v),1.1*max(theta*std_v)]);
        title(['Iteration',int2str(iter),':  reconstructed signal']);
        subplot(3,1,2), stem(indi,'.'); axis([0,M+1,-0.1,1.1]); title('Indicator of nonzero')
        subplot(3,1,3), stem(pi_tilde,'.'); axis([0,M+1,-0.1,1.1]);title('Probability of nonzero')
        xlabel('Coefficient Index'); 
        pause(0.001);
    end

end

m_theta=mean(cat(2,samples(:).theta),2);
