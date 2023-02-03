function [white_noise] = reshape_white_noise(white_noise,s,b,se)
% Reshape the white noise to the shape the algorithm need
% white_noise.all      1D white noise
% white_noise.a        [se.n_pts(1),se.n_pts(2)] white noise (same numbers
% as first entries of white_noise_all
% white_noise.b        [se.n_pts(1),se.n_pts(2)] white noise
% white_noise.c        [1] white noise
% white_noise.d        [b.n] white noise
% by Sebastian Reuschen

    white_noise.a = reshape(white_noise.all(1:se.n_pts(1)*se.n_pts(2)),[se.n_pts(1),se.n_pts(2)]);
    used_index = se.n_pts(1)*se.n_pts(2);

    if s.flag_kit == false
        white_noise.b = reshape(white_noise.all(used_index+1:used_index+se.n_pts(1)*se.n_pts(2)),[se.n_pts(1),se.n_pts(2)]);
        used_index = used_index+se.n_pts(1)*se.n_pts(2);
    end

    if s.flag_kit == 2
        white_noise.c = white_noise.all(used_index+1);
        used_index = used_index+1;
    end 

    if strcmp(b.model,'uncertain')
        white_noise.d(:,1) = white_noise.all(used_index+1:used_index+b.n);
        used_index = used_index+b.n;
    end

    if used_index ~= size(white_noise.all)
        throw(MException('ReshapeFUNC:BadRreshape','ERROR: The number of used random variables is not equal to the number of produced random variables.'))
    end
    


end