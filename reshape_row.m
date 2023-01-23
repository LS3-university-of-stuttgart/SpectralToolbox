function a = reshape_row(a)
% RESHAPE_ROW is a functional form of the COLON operator
% it returns any matrix reshaped to a row vector on-the-fly
% version 06 Mar 2008 / WN

n = numel(a);
a = reshape(a,1,n);