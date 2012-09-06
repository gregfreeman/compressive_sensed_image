function [m_theta, samples,burnin_samples]=TSWCSmcmc_NoSc(Phi, v, Ms, IdxParent, IdxChildren, wLevel,hyperpara, MCMCpara, plotflag)
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
PhiTPhi=sum(Phi.*Phi,1)';


% --------------
% Gibbs Sampling
% --------------
samples(MCMCpara.nCollect).theta=[];
burnin_samples=cell(MCMCpara.nBurnin,1);
for iter=1:(MCMCpara.nBurnin+MCMCpara.nCollect)

    % (1) theta -- sequentially drawn
    % \tilde{alpha}_{i}, M x 1
    alpha_tilde=alpha(wLevel)'+alpha0*PhiTPhi;
    % i=1
    A=Phi*theta;
    v_tilde=v-A+Phi(:,1)*theta(1);
    mu_tilde=(alpha0*Phi(:,1)'*v_tilde)/alpha_tilde(1);
    ratio=sqrt(alpha(wLevel(1))/alpha_tilde(1))*exp(0.5*alpha_tilde(1)*mu_tilde*mu_tilde)*PI(1)/(1-PI(1));
    theta_old=theta(1);
    if isinf(ratio)
        pi_tilde(1)=1;
        theta(1)=normrnd(mu_tilde, sqrt(1/alpha_tilde(1)));
    else
        pi_tilde(1)=ratio/(ratio+1);
        if rand<pi_tilde(1)
            theta(1)=normrnd(mu_tilde, sqrt(1/alpha_tilde(1)));
        else
            theta(1)=0;
        end
    end
    for i=2:M
        % update A
        A=A+Phi(:,i-1)*(theta(i-1)-theta_old);
        % \tilde{mu}_{i}, M x 1
        v_tilde=v-A+Phi(:,i)*theta(i);
        mu_tilde=(alpha0*Phi(:,i)'*v_tilde)/alpha_tilde(i);
        % \tilde{pi}_{i}, M x 1
        ratio=sqrt(alpha(wLevel(i))/alpha_tilde(i))*exp(0.5*alpha_tilde(i)*mu_tilde*mu_tilde)*PI(i)/(1-PI(i));
        % sample theta
        theta_old=theta(i);
        if isinf(ratio)
            pi_tilde(i)=1;
            theta(i)=normrnd(mu_tilde, sqrt(1/alpha_tilde(i)));
        else
            pi_tilde(i)=ratio/(ratio+1);
            if rand<pi_tilde(i)
                theta(i)=normrnd(mu_tilde, sqrt(1/alpha_tilde(i)));
            else
                theta(i)=0;
            end
        end
    end
    z=find(theta~=0);
    
    % (2) alpha
    for s=1:L
        idx=find(theta~=0 & wLevel==s);
        if ~isempty(idx)
            alpha(s)=gamrnd(c+0.5*length(idx), 1/(d+0.5*sum(theta(idx).^2,1)));
        end
    end
        
    % (3) PI
    indi=zeros(M,1); indi(z)=1;     % nonzero coefficient indicator
    Trans=zeros(2,2,L);             % store current transition matrix
    % for root level
    idx=find(wLevel==1); %idxLevelStart(1):idxLevelEnd(1);
    Nnon0=sum(indi(idx),1);
    beta_sample=betarnd(er+Nnon0,fr+length(idx)-Nnon0);
    PI(idx)=beta_sample;
    Trans(2,:,1)=[1-beta_sample, beta_sample];
    % for other levels
    for s=2:L
        idx=find(wLevel==s);%idxLevelStart(s):idxLevelEnd(s);
        % for parent=0
        idx0=indi(IdxParent(idx))==0;  % the indexes of 0 parent value in idx
        Nnon0=sum(indi(idx(idx0)),1);
        beta_sample=betarnd(e0(s)+Nnon0, f0(s)+length(idx0)-Nnon0);
        PI(idx(idx0))=beta_sample;
        Trans(1,:,s)=[1-beta_sample, beta_sample];
        % for parent~=0
        idx1=indi(IdxParent(idx))==1;  % the indexes of nonzero parent value in idx
        Nnon0=sum(indi(idx(idx1)),1);
        beta_sample=betarnd(e1(s)+Nnon0, f1(s)+length(idx1)-Nnon0);
        PI(idx(idx1))=beta_sample;
        Trans(2,:,s)=[1-beta_sample, beta_sample];
    end
    PI(PI==1)=1-eps;
    PI(PI==0)=eps;
    
    % (4) alpha0
    res=v-Phi(:,z)*theta(z);
    alpha0=gamrnd(a+N/2, 1/(b+res'*res/2));
    
    % Collect samples
    if iter<=MCMCpara.nBurnin
        burnin_samples{iter}= theta*std_v;
    end
    if iter>MCMCpara.nBurnin & (mod(iter-MCMCpara.nBurnin,MCMCpara.thinFactor)==1 | MCMCpara.thinFactor==1)
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
