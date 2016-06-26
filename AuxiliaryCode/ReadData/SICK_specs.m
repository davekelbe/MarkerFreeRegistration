function [ sampling, beam_diameter, min_obj_size, min_pts ] = SICK_specs( range, box_height)
%SICK_SPECS Computes the beam diameter and min_obj_size at a given range 
%   
%   Example: 
%   [sampling,beam_diameter,min_obj_size,min_pts] = SICK_specs(10,2);
%   
%   The beam diameter at 10m is .158m 
%   The minimum object size at 10m is 0.2016 diameter 
%
%
%   (C) David Kelbe, Rochester Institute of Technology 
    


beam_dia_at_optic = .008; % m
beam_divergence = .015; %radians 
angular_resolution = 0.25; % degrees

sampling = 2*tand(angular_resolution/2); %m

beam_diameter = [];
min_obj_size = [];
min_pts = inf;

if nargin<1;
    return
end
if isempty(range)
    return
end

if nargout > 1;
    beam_diameter = range * beam_divergence + beam_dia_at_optic;
end

if nargout > 2
    point_spacing = range*sampling;
    min_obj_size = point_spacing + beam_diameter;
end

if nargout >3
    min_pts = (box_height*min_obj_size)./point_spacing.^2;
end

end

