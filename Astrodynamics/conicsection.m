function v = conicsection(a, e, i, RAAN, w, step)
%% Compute conic section for orbital elements
%
%
%


%%
D2R = pi/180;
i = i*D2R;
RAAN = RAAN*D2R;
w = w*D2R;


%%
theta = [ 0:step:2*pi 0 ];
r = a*(1-e^2) ./ (1+e*cos(theta));
[x,y] = pol2cart(theta,r);

z = zeros(size(x));


%%
ECI_PQW = [cos(RAAN)*cos(w)-sin(RAAN)*sin(w)*cos(i) -cos(RAAN)*sin(w)-sin(RAAN)*cos(w)*cos(i) sin(RAAN)*sin(i);
           sin(RAAN)*cos(w)+cos(RAAN)*sin(w)*cos(i) -sin(RAAN)*sin(w)+cos(RAAN)*cos(w)*cos(i) -cos(RAAN)*sin(i);
           sin(w)*sin(i)                            cos(w)*sin(i)                             cos(i)];

       
%%
v = ECI_PQW*[x;y;z];