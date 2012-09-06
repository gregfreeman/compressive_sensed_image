%FastCSOperatorWrapper
% to use the FastCSOperator with the PDCO l1 optimization (slower) for
% direct comparison
% 
function y = FastCSOperatorWrapper(mode,m,n,x)

y = FastCSOperator(mode,m,n,x,true(1,n),n);
