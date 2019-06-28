function dx=nandiff(x)
% NANDIFF	Difference ignoring nans
% Like DIFF, but ignores NaNs. Instead of reporting NaN
% when difference includes a NaN, the differnence to next non-NaN number
% is logged. For matrices differences are calculated along columns. Array
% input or or input of dimension is not possible (yet).
%
% dx = nandiff(x)
%
% x	= vector or matrix input (MxN)
% dx	= vector or matrix of column differences ((M-1)xN)
%
% See also DIFF



[M,N]=size(x);	
dx=nan(M-1,N);

for j=1:N
  i=find(~isnan(x(:,j)));	% indices of nonan numbers in this column
  dx(i(1:end-1),j)=diff(x(i,j));
end

