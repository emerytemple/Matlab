function [x1,y1,z1, x2, y2, z2] = helix(LoD, D, theta, numpts, dir)
% LoD = L/D
% theta is in degrees
% dir = 1 for CCW, -1 for CW

	t = 0;
	theta = theta*(pi/180.0);
	dt = theta/(numpts-1);

	for i = 1:numpts
		% swap the sign inside both the sin() terms to rotate the other way
		x1(i) = D*cos(t);
		y1(i) = D*sin(dir*t);
		z1(i) = LoD*D*(t/theta);

		x2(i) = -D*cos(t);
		y2(i) = D*sin(-dir*t);
		z2(i) = LoD*D*(t/theta);		

		t = t + dt;
	end

	plot3(x1,y1,z1, x2, y2, z2)
end