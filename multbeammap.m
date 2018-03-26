function [MaxDif]=multbeammap(GEO_1,GEO_2)

% author: Seth Majeski
% Function: Determines The ground beam map based on the positions of Auris
% and a geostationary satellite
filename = '2week.csv';
numentries = 2000;            % number of auris positions
test1 = csvread(filename,1,1); 
test1(:,3)=test1(:,3)*1000;  %position of auris in lla (meters)
Earth = wgs84Ellipsoid();    %wgs84 model of earth
syms x y z                   %variables for finding position on earth
i=1;

P=zeros(numentries,3);
C=zeros(numentries,3);
llaPosition=zeros(numentries,3);
llaPosition1 = zeros(numentries,3);
p3=zeros(numentries,3);
p33=zeros(numentries,3);
 
llaGeosat = GEO_1;
llaGeosat2 = GEO_2;
%convert lla to ecef (m)
p11 = lla2ecef(llaGeosat2,'WGS84'); % geosat position in ecef
p1 = lla2ecef(llaGeosat, 'WGS84'); %  geosat position in ecef
p2 = lla2ecef(test1,'WGS84');    % example auris position in ecef
for i=1:numentries
p3(i,:) = p2(i,:)-p1; %difference in auris and geosat pos
p33(i,:)=p2(i,:)-p11;
end
%find the position by solving equations of derived line and earth ellipsoid
for i=1:numentries
% Solves the following equations for X, Y, Z

%                     LINE THROUGH GEOSAT/AURIS
% (X-geosat(x))/diff(x) = (Y-geosat(y))/diff(y) = (Z-geosat(z))/diff(z),
% where diff is the difference in Auris and geosat positions


%                      ELLIPSOID MODEL OF EARTH
% (X^2 + Y^2)/semimajoraxis^2 + Z^2/semiminoraxis^2 = 1


A = vpasolve((x-p1(1))/p3(i,1)==(y-p1(2))/p3(i,2),(x-p1(1))/p3(i,1)==(z-p1(3))/p3(i,3),(y-p1(2))/p3(i,2)==(z-p1(3))/p3(i,3),(x^2+ y^2)/Earth.SemimajorAxis^2 + z^2/Earth.SemiminorAxis^2 == 1); 
A1 = vpasolve((x-p11(1))/p33(i,1)==(y-p11(2))/p33(i,2),(x-p11(1))/p33(i,1)==(z-p11(3))/p33(i,3),(y-p11(2))/p33(i,2)==(z-p11(3))/p33(i,3),(x^2+ y^2)/Earth.SemimajorAxis^2 + z^2/Earth.SemiminorAxis^2 == 1);
%finds correct solution by checking which point is closer to the satellite
if(((A.x(1,1)-p2(i,1))^2+(A.y(1,1)-p2(i,2))^2+(A.z(1,1)-p2(i,3))^2) <((A.x(2,1)-p2(i,1))^2+(A.y(2,1)-p2(i,2))^2+(A.z(2,1)-p2(i,3))^2))
P(i,:) = [double(A.x(1,1)) double(A.y(1,1)) double(A.z(1,1))];
else
    P(i,:) = [double(A.x(2,1)) double(A.y(2,1)) double(A.z(2,1))];
end

if(((A1.x(1,1)-p2(i,1))^2+(A1.y(1,1)-p2(i,2))^2+(A1.z(1,1)-p2(i,3))^2) <((A1.x(2,1)-p2(i,1))^2+(A1.y(2,1)-p2(i,2))^2+(A1.z(2,1)-p2(i,3))^2))
C(i,:) = [double(A1.x(1,1)) double(A1.y(1,1)) double(A1.z(1,1))];
else
    C(i,:) = [double(A1.x(2,1)) double(A1.y(2,1)) double(A1.z(2,1))];
end

%convert to LLA
if(isreal(A.x(1,1))&&isreal(A.x(2,1)))
    llaPosition(i,:) = ecef2lla(P(i,:),'WGS84');
    llaPosition1(i,:) = ecef2lla(C(i,:),'WGS84');
end

end
%{
ecePos = lla2ecef(llaPosition);
[t1,t2,t3] = ellipsoid(0,0,0,Earth.SemimajorAxis,Earth.SemimajorAxis,Earth.SemiminorAxis);
figure
hold

%map globe
surf(t1,t2,t3);
scatter3(p1(1),p1(2),p1(3),'LineWidth',1.5);
plot3(p2(:,1),p2(:,2),p2(:,3),'LineWidth',1.5);
plot3(ecePos(:,1), ecePos(:,2), ecePos(:,3), 'black', 'LineWidth',2);
legend('','Geo-Satellite','Auris','Position on Earth')
%}
%map usa
DifDist = sqrt((C(1,:)-P(1,:)).^2 + (C(2,:)-P(2,:)).^2 + (C(3,:)-P(3,:)).^2);
MaxDif = max(DifDist);
fprintf('The maximum difference in the ground beam map positions: %11.2f m\n',MaxDif); 
figure
ax = usamap('conus');
oceanColor = [.5 .7 .9];
setm(ax, 'FFaceColor', oceanColor)
states = geoshape(shaperead('usastatehi', 'UseGeoCoords', true));
geoshow(states)
geoshow(llaPosition(:,1),llaPosition(:,2),'DisplayType','Point','Marker','*', 'MarkerSize',3)
geoshow(llaPosition1(:,1),llaPosition1(:,2),'DisplayType','Point','Marker','+','MarkerEdgeColor','blue', 'MarkerSize',3)
end