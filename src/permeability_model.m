function sample = permeability_model(white_noise_from_MCMC,parameters)
    white_noise.all = white_noise_from_MCMC;
    %Backtransformation of white noise
    white_noise = reshape_white_noise(white_noise,parameters.s,parameters.b,parameters.se);

    sample  = generate_randomfield(parameters.s,parameters.b,parameters.options,white_noise,parameters.FFTQe);
    sample = reshape(sample,[prod(parameters.s.n_pts),1]);
end