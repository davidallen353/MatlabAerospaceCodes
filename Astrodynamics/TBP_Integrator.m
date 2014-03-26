function [r,v,tout]=TBP_Integrator(r_0,v_0,tspan,mu)
% The purpose of this program is to, given an initial position and
% velocity, integrate the two body
%
% Inputs:
%   r_0: Initial position
%   v_0: Initial veloctity
%   tspan: The time span of integration. Must be a vector. If two elements,
%           the time steps will be chosed by the integrator. If more than
%           two elements, will report results at each time step
%   mu: Standard gravitational parameter. Earth if not specified
%
% Outputs:
%   r: Matrix of position vectors
%   v: Matrix of velocity vectors
%   tout: The time steps of each velocity vector
%
% Written by:
%   David Allen
%   davidallen@vt.edu
%   February 6, 2013

%========== Error Checking and data formating==============================
if nargin==3
    mu=3.986e5;
elseif nargin<3
    error('Too few input arguements')
elseif nargin>4
    error('Too many input arguements')
end

if min(size(r_0))~=1
    error('r_0 must be a vector')
elseif min(size(v_0))~=1
    error('v_0 must be a vector')
elseif min(size(tspan))~=1
    error('tspan must be a vector. If two elements, the time steps will be chosed by the integrator. If more than two elements, will report results at each time step')
elseif length(r_0)~=3
    error('r_0 must be a vector in 3-space')
elseif length(v_0)~=3
    error('v_0 must be a vector in 3-space')
elseif length(tspan)==1
    error('tspan must be a vector. If two elements, the time steps will be chosed by the integrator. If more than two elements, will report results at each time step')
end

r_0=[r_0(1);r_0(2);r_0(3)];
v_0=[v_0(1);v_0(2);v_0(3)];

x_0=[r_0;v_0];%stack the state variables
%==========================================================================

%======= Function to integrate=============================================
    function dx=TBP(t,x,mu)
        dx=zeros(6,1);
        
        dx(1)=x(4);
        dx(2)=x(5);
        dx(3)=x(6);
        
        R=norm([x(1),x(2),x(3)]);
        
        dx(4)=-(mu/R^3)*x(1);
        dx(5)=-(mu/R^3)*x(2);
        dx(6)=-(mu/R^3)*x(3);
    end
%==========================================================================

[tout,x]=ode45(@TBP,tspan,x_0,[],mu);

r=x(:,1:3)';
v=x(:,4:6)';

end
