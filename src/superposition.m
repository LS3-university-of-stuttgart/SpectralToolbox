function  [v,ve] = superposition(Qse1,u,H,n_pts,n_pts_e,nd,options,fftnQse1)
% SUPERPOSITION evaluates a superposition by user-specified method
% version 01 august 2007 / WN

if nargin < 8
  fftnQse1 = [];
end

switch options.superpos
  case 'standard'
    m         = numel(H);
    v         = zeros([n_pts 1]);
    for i=1:m
      v_i     = shiftaround(Qse1,H(i),n_pts,nd);
      v       = v + v_i*u(i);
    end
    ve=v;
  case 'fft'
    ui   = injection(u,H,n_pts);              % injection
    ue   = padding(ui,n_pts,n_pts_e,nd);      % embedding
    % convolution via FT
    ve   = fftn(ue);
    if ~isempty(fftnQse1)
      ve = ve .* fftnQse1;
    else
      ve = ve.*fftn(Qse1);
    end
    ve   = ifftn(ve);
    ve   = real(ve);
    v    = extraction(ve,n_pts,nd);           % extraction
  otherwise
    error('SUPERPOSITION.input.options: superposition option must be STANDARD or FFT')
end