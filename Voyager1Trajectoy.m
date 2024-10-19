% This is the VEEGA transfer scenario (2323 days)

%Edited for Voyager1

clear all
clc
close all

% Start Data @Earth for Sept 05 1977
% source: https://ssd.jpl.nasa.gov/?horizons
%SR
StartPose = 1.495978707e11.*[ 9.646121667104959E-01 -3.052734146049315E-01 -7.160201343849379E-05]; %1.495978707e11.* conversion from AU to meters. infor from X in $$SOE
StartVel = [4.868314521712150E-03 1.635386972936581E-02 1.009178167790655E-06];   %AU/day Info from VX in $$SOE

tau1 = 546*86400; %days*86400 sec SR
exitAng = 71*pi/180;
exitVec = 10240.*[cos(exitAng) sin(exitAng) 0];
%JENTRY = norm(exitVec);%q2


r0 = StartPose;
r0(3) = 0;  % flatten the solar system for conveniance
v0 = 1.73145684e6.*[StartVel(1) StartVel(2) 0] + exitVec; %AU per day in m/s multiply by 1.73145684e6
[rx0, ry0, rz0,vx0,vy0,vz0,tend] = RKballistic3D(tau1,v0,r0);
JENTRY=sqrt(vx0(end)^2+vy0(end)^2);
%jENTRYx=norm(v0);
%fprintf('\nHeliocentric velocity after earth is %s m/s\n\n', num2str(norm(v0))); %Q2

figure(1)
plot(rx0,ry0,'r', 'LineWidth',1)
hold on

%% Jupiter Encounter March 5th, 1979
JEnc = 1.495978707e11.*[-3.213666304496740E+00 4.193866288474054E+00 5.469815440880738E-02];
JEncVel = [-6.074779618322786E-03 -4.242598233389083E-03 1.535285223164491E-04];

JR = sqrt(dot(JEnc,JEnc));

tau2 = 618*86400; %days*86400 sec SR
exitAng = 197.195*pi/180;
exitVec = 10975.*[cos(exitAng) sin(exitAng) 0];
JEXIT = exitVec;%q2
JV=JEncVel;%q2

r0 = JEnc;
r0(3) = 0;  % flatten the solar system for conveniance
v0 = 1.73145684e6.*[JEncVel(1) JEncVel(2) 0] + exitVec; %AU per day in m/s multiply by 1.73145684e6
[rx1, ry1, rz1, tend] = RKballistic3D(tau2,v0,r0);
jEXITx=norm(v0);

%fprintf('\nHeliocentric velocity after Jupiter is %s m/s\n\n', num2str(norm(v0))); %Q2

figure(1)
plot(rx1,ry1,'r', 'LineWidth',1)
hold on

%% Saturn Encounter  November 12, 1980
SEnc = 1.495978707e11.*[-9.492363608272317E+00 -3.575002051307236E-01 3.834817000651295E-01];
SEncVel = [-8.852696373106042E-05 -5.585528597312641E-03 1.014779242821678E-04];

tau3 = 731*86400;
exitAng = 35*pi/180;
exitVec = -12050.*[cos(exitAng) sin(exitAng) 0];

r0 = SEnc;
r0(3) = 0;  % flatten the solar system for convenience
v0 = 1.73145684e6.*[SEncVel(1) SEncVel(2) 0] + exitVec; %AU per day in m/s multiply by 1.73145684e6
[rx2, ry2, rz2, tend] = RKballistic3D(tau3,v0,r0);

figure(1)
plot(rx2,ry2,'r', 'LineWidth',1)
hold on


SR = sqrt(dot(StartPose,StartPose));

% End Data @Jupiter for April 04, 2028
EndPose = 1.495978707e11.*[-5.403570455360613E+00 5.759180903371619E-01 1.185602847664122E-01];
EndVel = [-8.908379419632020E-04 -7.152519660681221E-03 4.963758207185361E-05];

TR = sqrt(dot(SEnc,SEnc));

%% Add planetary orbits
angle = linspace(0,2*pi,361);
for i=1:length(angle)
    ERorbit(i,1) = SR.*cos(angle(i));
    ERorbit(i,2) = SR.*sin(angle(i));
    ERorbit(i,3) = StartPose(3);
    
    Jorbit(i,1) = JR.*cos(angle(i));
    Jorbit(i,2) = JR.*sin(angle(i));
    Jorbit(i,3) = JEnc(3);
    
    Sorbit(i,1) = TR.*cos(angle(i));
    Sorbit(i,2) = TR.*sin(angle(i));
    %Sorbit(i,3) = EndPose(3);
end


for i = 1:length(rx0)
    OrbDist(i,1) = sqrt(rx0(i)^2+ry0(i)^2); %distance at each point in the orbit
    OrbDist(i,2) = (i-1)*1000;
    Exp(i) = 1000./((OrbDist(i,1)/1.495978707e11)^2);
end

offset = length(rx0);

for i = 1:length(rx1)
    OrbDist(i+offset,1) = sqrt(rx1(i)^2+ry1(i)^2); %distance at each point in the orbit
    OrbDist(i+offset,2) = (i-1)*1000+tau1;
    Exp(i+offset) = 1000./((OrbDist(i+offset,1)/1.495978707e11)^2);
end

offset = offset+length(rx1);

for i = 1:length(rx2)
    OrbDist(i+offset,1) = sqrt(rx2(i)^2+ry2(i)^2); %distance at each point in the orbit
    OrbDist(i+offset,2) = (i-1)*1000+tau1+tau2;
    Exp(i+offset) = 1000./((OrbDist(i+offset,1)/1.495978707e11)^2);
end

offset = offset+length(rx2);

TotalExpVOYAGER = sum(Exp)./86400;   %Total 1-AU days of exposure in the DIRECT Trajectory
figure(1)
%plot(ERorbit(:,1),ERorbit(:,2),'b:',StartPose(1),StartPose(2),'go',Sorbit(:,1),Sorbit(:,2),'b:',Jorbit(:,1),Jorbit(:,2),'b:',JEnc(1),JEnc(2),'ro',SEnc(1),SEnc(2),'wo', 'LineWidth',2)
plot(ERorbit(:,1),ERorbit(:,2),'b:',StartPose(1),StartPose(2),'co',Sorbit(:,1),Sorbit(:,2),'b:',Jorbit(:,1),Jorbit(:,2),'b:',SEnc(1),SEnc(2),'k_', 'LineWidth',2);
title('Reconstructed trajectory of Voyager 1')
plot(JEnc(1),JEnc(2),'r*', 'LineWidth',2);
plot(JEnc(1),JEnc(2),'o', 'Color', [1 0.5 0], 'LineWidth',3);
plot(SEnc(1),SEnc(2),'yo', 'LineWidth',1)
plot(SEnc(1),SEnc(2),'w_', 'LineWidth',1)
plot(StartPose(1),StartPose(2),'g*', 'LineWidth',2)

set(gca,'color',[0 0 .2])
hold on
plot(0, 0, 'hexagram', 'MarkerSize', 10, 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'none'); %Mr.Sun, Sun, Mr.Golden Sun

hold off
CFig = get(1);
chil = CFig.Children(1);
set(chil,'FontSize',16);
set(chil,'LineWidth',2);

figure(2)
plot(OrbDist(:,2)./(86400*365.24),OrbDist(:,1)./1.496e11,'LineWidth',2)
set(gca,'color',[0 0 0])
axis tight
title('Heliocentric Distance versus Time','FontSize',24)
ylabel('Distance in AU','FontSize',18)
xlabel('Time in Years','FontSize',18)
CFig = get(2);
chil = CFig.Children(1);
set(chil,'FontSize',16);
set(chil,'LineWidth',2);

disp('~~ Orbit Reconstruction Complete ~~')