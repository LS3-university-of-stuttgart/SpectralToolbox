function [e_i] = unit_vector(n,i)
% UNIT_VECTOR returns all-zero vector with unit entry at position i
% version 01 august 2007 / WN

e_i    = zeros([n,1]);
e_i(i) = 1;
