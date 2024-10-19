function [rx, ry, rz, tend] = RKballistic3D(tau,v0,r0)

%% CONSTANTS

G = 6.67e-11;       % gravitational constant
M = 1.989e30;       % mass of the Sun
dt = 1000;

%% INITIAL CONDITIONS
vx(1) = v0(1);  % initial velocity
vy(1) = v0(2);
vz(1) = v0(3);
vm(1) = sqrt(vx(1).^2+vy(1).^2+vz(1).^2);

rx(1) = r0(1);  %initial position
ry(1) = r0(2);
rz(1) = r0(3);
rm(1) = sqrt(rx(1).^2+ry(1).^2+rz(1).^2);  %magnitude of the position vector

am(1) = -1*G*M./rm(1).^2;     % magnitude of acceleration  vector
ax(1) = am(1).*rx(1)./rm(1);
ay(1) = am(1).*ry(1)./rm(1);
az(1) = am(1).*rz(1)./rm(1);
t(1) = 0;

%% Advance the Payload Position

i= 1;

while(t < tau)
    % Update Time
    i = i+1;
    t(i) = (i-1)*dt;
    tend = t(i);    % output value
    % Interleaved RK4 Method
    k1rx = vx(i-1);
    k1ry = vy(i-1);
    k1rz = vz(i-1);
    j1vx = -1*G*M*rx(i-1)/(rm(i-1)^3);
    j1vy = -1*G*M*ry(i-1)/(rm(i-1)^3);
    j1vz = -1*G*M*rz(i-1)/(rm(i-1)^3);
    
    k2rx = vx(i-1)+0.5*j1vx;
    k2ry = vy(i-1)+0.5*j1vy;
    k2rz = vz(i-1)+0.5*j1vz;
    j2vx = -1*G*M*(rx(i-1)+0.5*k1rx)/(sqrt((rx(i-1)+0.5*k1rx).^2+(ry(i-1)+0.5*k1ry).^2+(rz(i-1)+0.5*k1rz).^2))^3;
    j2vy = -1*G*M*(ry(i-1)+0.5*k1ry)/(sqrt((rx(i-1)+0.5*k1rx).^2+(ry(i-1)+0.5*k1ry).^2+(rz(i-1)+0.5*k1rz).^2))^3;
    j2vz = -1*G*M*(rz(i-1)+0.5*k1rz)/(sqrt((rx(i-1)+0.5*k1rx).^2+(ry(i-1)+0.5*k1ry).^2+(rz(i-1)+0.5*k1rz).^2))^3;
    
    k3rx = vx(i-1)+0.5*j2vx;
    k3ry = vy(i-1)+0.5*j2vy;
    k3rz = vz(i-1)+0.5*j2vz;
    j3vx = -1*G*M*(rx(i-1)+0.5*k2rx)/(sqrt((rx(i-1)+0.5*k2rx).^2+(ry(i-1)+0.5*k2ry).^2+(rz(i-1)+0.5*k2rz).^2))^3;
    j3vy = -1*G*M*(ry(i-1)+0.5*k2ry)/(sqrt((rx(i-1)+0.5*k2rx).^2+(ry(i-1)+0.5*k2ry).^2+(rz(i-1)+0.5*k2rz).^2))^3;
    j3vz = -1*G*M*(rz(i-1)+0.5*k2rz)/(sqrt((rx(i-1)+0.5*k2rx).^2+(ry(i-1)+0.5*k2ry).^2+(rz(i-1)+0.5*k2rz).^2))^3;
    
    k4rx = vx(i-1)+j2vx;
    k4ry = vy(i-1)+j2vy;
    k4rz = vz(i-1)+j2vz;
    j4vx = -1*G*M*(rx(i-1)+k3rx)/(sqrt((rx(i-1)+k3rx).^2+(ry(i-1)+k3ry).^2+(rz(i-1)+k3rz).^2))^3;
    j4vy = -1*G*M*(ry(i-1)+k3ry)/(sqrt((rx(i-1)+k3rx).^2+(ry(i-1)+k3ry).^2+(rz(i-1)+k3rz).^2))^3;
    j4vz = -1*G*M*(rz(i-1)+k3rz)/(sqrt((rx(i-1)+k3rx).^2+(ry(i-1)+k3ry).^2+(rz(i-1)+k3rz).^2))^3;
    
    rx(i) = rx(i-1)+(dt/6)*(k1rx+2*k2rx+2*k3rx+k4rx);
    ry(i) = ry(i-1)+(dt/6)*(k1ry+2*k2ry+2*k3ry+k4ry);
    rz(i) = rz(i-1)+(dt/6)*(k1rz+2*k2rz+2*k3rz+k4rz);
    vx(i) = vx(i-1)+(dt/6)*(j1vx+2*j2vx+2*j3vx+j4vx);
    vy(i) = vy(i-1)+(dt/6)*(j1vy+2*j2vy+2*j3vy+j4vy);
    vz(i) = vz(i-1)+(dt/6)*(j1vz+2*j2vz+2*j3vz+j4vz);
    
    % Calculate Magnitudes of position, velocity and acceleration vectors
    rm(i) = sqrt(rx(i).^2+ry(i).^2);  
    vm(i) = sqrt(vx(i).^2+vy(i).^2);          
    am(i) = -1*G*M./rm(i).^2;        
    ax(i) = am(i).*rx(i)./rm(i);
    ay(i) = am(i).*ry(i)./rm(i);
    az(i) = am(i).*rz(i)./rm(i);
end