function rotated_coords = rotation(input_XY,center,anti_clockwise_angle)

[r,c] = size(center);

% Format the coordinate of the center of rotation
center_coord = input_XY;
center_coord(:,1) = center(1);
center_coord(:,2) = center(2);

% Turns the angles given to be such that the +ve is anti-clockwise and -ve is clockwise
anti_clockwise_angle = -1*anti_clockwise_angle;
% if in degrees, convert to radians because that's what the built-in functions use. 

anti_clockwise_angle = deg2rad(anti_clockwise_angle);


%Produce the roation matrix
rotation_matrix = [cos(anti_clockwise_angle),-1*sin(anti_clockwise_angle);...
                   sin(anti_clockwise_angle),cos(anti_clockwise_angle)];
%Calculate the final coordinates
rotated_coords = ((input_XY-center_coord) * rotation_matrix)+center_coord;

