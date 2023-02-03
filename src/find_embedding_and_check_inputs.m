function [se] = find_embedding_and_check_inputs(s,options,b)
% FIND_EMBEDDING chooses minimum embedding size given grid and geostatistical model
% by Sebastian Reuschen
    % checking required input parameters
    s = check_structure_s(s);
    s = ndgrid_setup(s);
    b = check_structure_b(b);

    se         = find_embedding(s,options);

end