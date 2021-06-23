function [Res] = isPositive(p1,p2,p3)
% To=Tindex;
% p1=Tria(To,1);
% p2=Tria(To,2);
% p3=Tria(To,3);

% xy_p1=Tria.Points(p1,:);
% %xy_p1=[xy_p1 0];
% xy_p2=Tria.Points(p2,:);
% %xy_p2=[xy_p2 0];
% xy_p3=Tria.Points(p3,:);
% %xy_p3=[xy_p3 0];

xy_p1=p1;
xy_p2=p2;
xy_p3=p3;

p1_p2=[xy_p2(1)-xy_p1(1),xy_p2(2)-xy_p1(2)];
p2_p3=[xy_p3(1)-xy_p2(1),xy_p3(2)-xy_p2(2)];
p3_p1=[xy_p1(1)-xy_p3(1),xy_p1(2)-xy_p3(2)];

v1=[p1_p2 0];
v2=[p2_p3 0];
v3=[p3_p1 0];

cross1=cross(v1,v2);
cross2=cross(v2,v3);
cross3=cross(v3,v1);

%O=[cross1(3) cross2(3) cross3(3)];

Res=true;
if cross1(3)<0
    Res=false;
end

end