function [n_white_noise] = compute_n_white_noise(se,s,b)
% Computes the number of white noise that needs to be generated for the
% model.
% by Sebastian Reuschen
    n_white_noise = prod(se.n_pts);
    if s.flag_kit == false
        n_white_noise = n_white_noise + prod(se.n_pts);
    end

    if s.flag_kit == 2
        n_white_noise = n_white_noise + 1;
        %white_noise(3,:,:)=randn(1,1);
    end 

    if strcmp(b.model,'uncertain')
        n_white_noise = n_white_noise + b.n;
        %white_noise(4,:,:)=randn(b.n,1);
    end


end