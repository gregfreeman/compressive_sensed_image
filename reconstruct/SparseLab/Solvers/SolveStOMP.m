function [sol, numIters,nActiveSet,normRes] = SolveStOMP(A, y, N, thresh, param, maxIters, verbose, OptTol)
% SolveStOMP: Implementation of Iterative Threshold-Selective Projection
% algorithm
% Usage
%	[sol, numIters] = SolveStOMP(A, y, N, thresh, param, maxIters,
%	verbose, OptTol)
% Input
%	A           Either an explicit nxN matrix, with rank(A) = min(N,n) 
%               by assumption, or a string containing the name of a 
%               function implementing an implicit matrix (see below for 
%               details on the format of the function).
%	y           vector of length n.
%   N           length of solution vector. 
%   thresh      Thresholding strategy: FDR or FAR. default is FDR.
%   param       Sensitivity parameter for threshold selection.
%	maxIters    maximum number of StOMP iterations to perform, default 10. 
%   verbose     1 to print out detailed progress at each iteration, 0 for
%               no output (default).
%	OptTol      Error tolerance, default 1e-5.
% Outputs
%	 sol        Solution of StOMP.
%	 numIters   Total number of steps taken.
% Description
%   SolveStOMP implements the Stagewise Ortogonal Matching Pursuit, 
%   as described in the paper
%   D.L. Donoho et al., "Sparse Solution of Underdetermined Linear 
%   Equations by Stagewise Ortogonal Matching Pursuit". 
%   It essentially computes an approximate solution to the problem 
%     min ||x||_1  s.t.  Ax = y
%   The matrix A can be either an explicit matrix, an implicit operator
%   implemented as an m-file, or a matlab object. If using the implicit form, 
%   the user should provide the name of a function of the following format:
%     y = OperatorName(mode, m, n, x, I, dim)
%   This function gets as input a vector x and an index set I, and returns
%     y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2. 
%   A is the m by dim implicit matrix implemented by the function. I is a
%   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
%   vector of length n is mode = 1, or a vector of length m is mode = 2.
%   If using the object form, the matlab object should have a member
%   function named csGenOperator with the following format:
%     y = csGenOperator(this, mode, m, n, x, I, dim)
%   Where this is an object and other parameters conform to functional mode
%   This function is used as input to Michael Saunders' least-squares
%   solver, LSQR, available at:
%     http://www.stanford.edu/group/SOL/software/lsqr.html
% See Also
%   SolveLasso, SolveBP, SolveOMP
%

dispErr = false;
cheat = false;

if nargin < 8,
    OptTol = 1e-5;
end
if nargin < 7,
    verbose = 0;
end
if nargin < 6,
    maxIters = 10;
end
if nargin < 4,
    thresh = 'FDR';
    param = 0.5;
end

% Initialize threshold parameters
switch lower(thresh)
    case 'fdr'
        q = param;
    case 'far'
        alpha_0 = param;
    case 'rcl'
        maxScale = param;
end

explicitA = ~(ischar(A) || isa(A, 'function_handle') || isobject(A));
objectA = isobject(A);
n = length(y);

if (dispErr && objectA),
   corrFig = figure();
end

% Set LSQR parameters
damp   = 0;
atol   = 1.0e-6;
btol   = 1.0e-6;
conlim = 1.0e+10;
itnlim = n;
show   = 0;

% Initialize
iter = 1;
res = y;  
normy = norm(y);
activeSet = [];
x_I = [];
Ifull = 1:N;
done = 0;

sol = zeros(N,1);

while ~done
    % Compute residual correlations
    if (explicitA)
        corr = sqrt(n) .* A'*res ./ norm(res);
    else
        if objectA,
             corr = csGenOperator(A,2,n,N,res,Ifull,N);
             %apply perceptual basis according to object parameters
             corr = adjustEstimate(A,corr);
        else
             corr = feval(A,2,n,N,res,Ifull,N);
        end
        
        corr = sqrt(n) .* corr ./ norm(res);
        
    end
    
    switch lower(thresh)
        case 'fdr'
            thr = fdrthresh(corr, q);
        case 'far'
            thr = norminv2(1 - alpha_0/2, 0, 1);
        case 'rcl'
            maxCorr = max(abs(corr));
            medCorr = median(abs(corr));
            thr = maxCorr * maxScale;
            if thr < (medCorr / maxScale)
                thr = maxCorr * 1.1;
            end
    end
    
    if (cheat && objectA),
        trueI = zeros(N,1);
        trueI(getActiveInds(A)) = 1;
        corrNoise = corr(trueI == 0);
        corrNoise = reverse(sort(abs(corrNoise)));
        thr = corrNoise(1)*0.7;
    end

    % Apply hard thresholding
    I = find(abs(corr) > thr);
    
    % If no new entries, we are done
    J = union(activeSet, I);
    
    if (dispErr && objectA),
        trueI = zeros(N,1);
        calcI = zeros(N,1);
        trueI(getActiveInds(A)) = 1;
        calcI(J) = 1;
        hits = sum(trueI.*calcI);
        falsePos = sum(calcI - trueI.*calcI);
        corrSig = corr(trueI == 1);
        corrSig = reverse(sort(abs(corrSig)));
        corrNoise = corr(trueI == 0);
        corrNoise = reverse(sort(abs(corrNoise)));
        fprintf('Iteration %d: Hits = %d, falsePos = %d\n', iter, hits, falsePos);
        figure(corrFig);
        plot(1:length(corrSig),corrSig,'b.',1:length(corrNoise),corrNoise,'r.');
    end
    if (length(J) == length(activeSet)) 
        done = 1;
    else
        activeSet = J;
        
        % Compute current estimate and residual
        if (explicitA)
            x_I = A(:, activeSet) \ y;
            res = y - A(:, activeSet)*x_I;
        else
            % Solve (A_I^T*A_I)x_I = y using lsqr
            p = length(activeSet);
            %[x_I, istop, itn] = lsqr_ms(n, p, A, activeSet, N, y, damp, atol, btol, conlim, itnlim, show);
            [x_I, flag] = lsqr(@lsqrMat,y,OptTol,20);

            % Compute residual
            if objectA,
                Ax_I = csGenOperator(A,1,n,p,x_I,activeSet,N);
            else
                Ax_I = feval(A,1,n,p,x_I,activeSet,N);
            end
            res = y - Ax_I;
        end

        if (verbose)
            fprintf('Iteration %d: |I| = %d, ||r||_2 = %g\n', iter, length(activeSet), norm(res));
            if (objectA)
                tempSol = zeros(N,1);
                tempSol(activeSet) = x_I;
                currErr = calcError(A,tempSol);
                fprintf('Current error: %g\n',currErr);
            end
        end
    end
    
    iter = iter+1;
    
    % Check stopping criteria
    if (iter > maxIters)
        done = 1;
    end
    if norm(res) <= OptTol*normy
        done = 1;
    end
    if (dispErr && objectA),
        sol(activeSet) = x_I;
        [A curErr] = visError(A,sol);
    end
end

sol = zeros(N,1);
sol(activeSet) = x_I;
numIters = iter;
nActiveSet=length(activeSet);
normRes=norm(res);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function y = lsqrMat(x,transpose)
        switch transpose
            case 'transp'
                if objectA,
                    y = csGenOperator(A,2,n,p,x,activeSet,N);
                else,
                    y = feval(A,2,n,p,x,activeSet,N);
                end
            case 'notransp'
                if objectA,
                    y = csGenOperator(A,1,n,p,x,activeSet,N);
                else,
                    y = feval(A,1,n,p,x,activeSet,N);
                end
        end
    end

end

%
% Copyright (c) 2006. Yaakov Tsaig
%  
            
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
