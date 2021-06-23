function [Tria,Point,Adja,Neig,Nmat,min,max,C,tou] = GenMesh(R,l,h)

%R=3; %circle radius
%l=2; %side of size l of the equilateral triangle

%Boundary of Omega 
%CIR=[R*cos([0:pi/100:2*pi]);R*sin([0:pi/100:2*pi])];
%h=4;%number of divided segments in a triangle edge

d=l/h;%distance between two consecutive poins in the triangle


min=3*(h-1)+3;%capturing the number of points in the trangle
%Edges of the triangle
e1=[[-l/2:d:l/2];(-l/(sqrt(3)*2))*ones(1,h+1)];
e2=[[0:(l/2)/h:l/2];-sqrt(3)*[0:(l/2)/h:l/2]+l/sqrt(3)];
e3=[[-l/2:(l/2)/h:0];sqrt(3)*[-l/2:(l/2)/h:0]+l/sqrt(3)];

N1=size(e1);
N2=size(e2);
N3=size(e3);

r=floor((R-(l/sqrt(3)))/d);%number of circles inside the circle of radius R
%circumference passing through the vertices of the triangle
CT=[l/sqrt(3)*cos([-5*pi/6:(2*pi/3)/(h):7*pi/6]);l/sqrt(3)*sin([-5*pi/6:(2*pi/3)/(h):7*pi/6])];
%circumference aproximation of Omega(radius R)
CRe=[R*cos([-5*pi/6:pi/(pi/acos((18*R*R-l*l)/(18*R*R))):7*pi/6]);R*sin([-5*pi/6:pi/(pi/acos((18*R*R-l*l)/(18*R*R))):7*pi/6])];

%%
[f,c]=size(CRe);
p0=[CRe(1,1) CRe(2,1)];
pn=[CRe(1,c) CRe(2,c)];
%%distance between the first and the last point in CRe
d2=sqrt((pn(1)-p0(1))^2+(pn(2)-p0(2))^2);

if 1 %R>2 || R<1
    if d2<=d/2%if the distance is <=d/2 the last point is not taken into account
        CRe=CRe(:,1:c-1);
    end
end    

%%
%hold on
plot(e1(1,:),e1(2,:),e2(1,:),e2(2,:),e3(1,:),e3(2,:));%,CT(1,:),CT(2,:));
a_x=[];
a_y=[];
R1=l/sqrt(3)+d;
%internal circumferences
for i=1:1:r
    %CR contains the points with radius R1
    CR=[R1*cos([-5*pi/6:pi/(pi/acos((18*R1*R1-l*l)/(18*R1*R1))):7*pi/6]);R1*sin([-5*pi/6:pi/(pi/acos((18*R1*R1-l*l)/(18*R1*R1))):7*pi/6])];
    [f,c]=size(CR);
    p0=[CR(1,1) CR(2,1)];
    pn=[CR(1,c) CR(2,c)];
    %distance between the first and the last point in CR
    d2=sqrt((pn(1)-p0(1))^2+(pn(2)-p0(2))^2);

    %display(R-R1)
    %If the last circumference is very close to the radius R, it is not taken into account
    if R-R1<0.1 && i==r %valor a la derecha de < esta puesto al azar
        break
    end
    if 1%R>2 || R<1
        %% if the distance is <=d/2 the last point in a circunference is not taken into account
        if d2<=d/2 %este d2<=dr no drbr cumplirse en 114
            a_x=[a_x CR(1,1:c-1)];
            a_y=[a_y CR(2,1:c-1)];
            %plot(CR(1,1:c-1),CR(2,1:c-1));
            R1=R1+d;
            continue
        end
    end
    %% if there is no problem all the points are created and added
    a_x=[a_x CR(1,:)];
    a_y=[a_y CR(2,:)];
    %plot(CR(1,:),CR(2,:));
    R1=R1+d;
end 
%%
a_x2=[e1(1,:) e2(1,N2(1,2)-1:-1:1) e3(1,N3(1,2)-1:-1:2) CT(1,2:h) CT(1,h+2:2*h) CT(1,2*h+2:3*h) a_x CRe(1,:)];
a_y2=[ e1(2,:) e2(2,N2(1,2)-1:-1:1) e3(2,N3(1,2)-1:-1:2) CT(2,2:h) CT(2,h+2:2*h) CT(2,2*h+2:3*h) a_y CRe(2,:)];

%fist point in the circumference with R
Nmat=length([CT(1,2:h) CT(1,h+2:2*h) CT(1,2*h+2:3*h) a_x]);
max=Nmat+min+1;


%%index boundary constraints
[~,cax]=size(a_x2);
[~,colcons]=size([e1(1,:) e2(1,N2(1,2)-1:-1:1) e3(1,N3(1,2)-1:-1:2)]);
[~,cpo]=size([e1(1,:) e2(1,N2(1,2)-1:-1:1) e3(1,N3(1,2)-1:-1:2) CT(1,2:h) CT(1,h+2:2*h) CT(1,2*h+2:3*h) a_x]);
innercons = [(1:colcons-1)' (2:colcons)'; colcons 1;];
outercons = [(cpo+1:cax-1)' (cpo+2:cax)'; cax cpo+1];
C = [innercons; outercons];
%First Triangulation with Constraints
DT=delaunayTriangulation(a_x2(:),a_y2(:),C);
%isInterior 
TF = isInterior(DT);
 
%%the final tringulation
Tria = triangulation(DT.ConnectivityList(TF,:),DT.Points(:,1),DT.Points(:,2));
triplot(Tria)
xlabel('x'), ylabel('y')
%triplot(DT.ConnectivityList(TF,:),DT.Points(:,1),DT.Points(:,2))  
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a_x=[e1(1,:) e2(1,N2(1,2)-1:-1:1) e3(1,N3(1,2)-1:-1:2) CT(1,2:h) CT(1,h+2:2*h) CT(1,2*h+2:3*h) a_x CRe(1,:)];
%a_y=[ e1(2,:) e2(2,N2(1,2)-1:-1:1) e3(2,N3(1,2)-1:-1:2) CT(2,2:h) CT(2,h+2:2*h) CT(2,2*h+2:3*h) a_y CRe(2,:)];



f2 = figure;
triplot(Tria)
xlabel('x'), ylabel('y')
IC = incenter(Tria);
Point=[a_x2(:) a_y2(:)];
 hold on
for i=1: size(Point,1)
    txt=['$$\bullet  ' int2str(i) '$$'];
    text(Point(i,1)-0.005,Point(i,2)-0.005,txt,'Interpreter','latex', ...
        'color','blue','FontSize',14);
    Adja(i).vect=[];
    for j=1:size(Tria,1)
        if sum(i==Tria(j,:))>0
            Adja(i).vect=[Adja(i).vect j];
        end
    end
end
count=0;
for i=1:size(Tria,1)
    txt=['$$\bullet  ' int2str(i) '$$'];
    text(IC(i,1),IC(i,2),txt,'Interpreter','latex','color','red','FontSize',14);  
    II=Tria(i,[1 2]);
    True=1;
    j=1;
    while True & j<(size(Tria,1)+1)
        if i~=j & (II==Tria(j,[1 2]) | II==Tria(j,[2 1]) | ...
            II==Tria(j,[2 3]) | II==Tria(j,[3 2]) | ...    
            II==Tria(j,[3 1]) | II==Tria(j,[1 3]) ) 
            Neig(i,1)=j;
            True=0;
        end
        j=j+1;
    end
    if True==1
        Neig(i,1)=-1;
    end
    II=Tria(i,[2 3]);
    True=1;
    j=1;
    while True & j<(size(Tria,1)+1) 
        if i~=j & (II==Tria(j,[1 2]) | II==Tria(j,[2 1]) | ...
            II==Tria(j,[2 3]) | II==Tria(j,[3 2]) | ...    
            II==Tria(j,[3 1]) | II==Tria(j,[1 3]) )
            Neig(i,2)=j;
            True=0;
        end
        j=j+1;
    end
    if True==1
        Neig(i,2)=-1;
    end   
    II=Tria(i,[3 1]);
    True=1;
    j=1;
    while True & j<(size(Tria,1)+1)
        if i~=j & (II==Tria(j,[1 2]) | II==Tria(j,[2 1]) | ...
            II==Tria(j,[2 3]) | II==Tria(j,[3 2]) | ...    
            II==Tria(j,[3 1]) | II==Tria(j,[1 3]) )
            Neig(i,3)=j;
            True=0;
        end
        j=j+1;
    end
    if True==1
        Neig(i,3)=-1;
    end
    A=Tria(i,1);
    B=Tria(i,2);
    C=Tria(i,3);
    Axy=Tria.Points(A,:);
    Bxy=Tria.Points(B,:);
    Cxy=Tria.Points(C,:);
    if isPositive(Axy,Bxy,Cxy)==true
        count=count+1;
    end
end
tou=false;
if(count==size(Tria,1))
    tou=true;
end  
end
% Label the vertices.

% hold on
% numvx = size(a_x,2);
% vxlabels = arrayfun(@(n) {sprintf('%d', n)}, (1:numvx)');
% Hpl = text(a_x(:)+0.05, a_y(:)+0.05, vxlabels, 'FontWeight', ...
%   'bold', 'HorizontalAlignment','center', 'BackgroundColor', ...
%   'none');
% hold off