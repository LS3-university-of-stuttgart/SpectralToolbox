function uQu = quadratic_form_FFT(FFTQe,u)
% QUADRATIC_FORM_FFT performs double convolution between Q and u via FFT
% to evaluate a quadratic matrix form u'*Q*u.
% Q represents the first row of a Toeplitz matrix, Qe the same
% embedded in a larger circulant matrix, and FFTQse is its FFT.
% U is the vector to be multiplied by Q.
% version 29 january 2008 / WN
%
% input parameters:
% FFTQe            FFT of the embedded version of Q
%                  must be larger than or equal in size to U.
%                  must be reshaped to the physical dimensionality.
% u                vector to be multiplied
%                  must be reshaped to the physical dimensionality.
%
% output parameters:
% uQu              result of quadratic form, a scalar

ue            = padding(u,size(FFTQe));
Ue            = fftn(ue);
uQu           = real(abs(Ue).^2.*FFTQe);
uQu           = sum(uQu(:))/numel(FFTQe);