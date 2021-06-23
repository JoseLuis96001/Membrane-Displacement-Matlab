function [Tria2]=refinement(Tria,C)

%[Tria,a_x2,a_y2,C,Nmat,min,max]=GenMesh2(0.8,0.7,7)
a_x2=Tria.Points(:,1)';
a_y2=Tria.Points(:,2)';

a_x=[];
a_y=[];
[n,~]=size(Tria.ConnectivityList);

for i=1:1:n
    [P]=Tria.ConnectivityList(i,:);
    P1=Tria.Points(P(1),:);
    P2=Tria.Points(P(2),:);
    P3=Tria.Points(P(3),:);
    mid1 = (P1(:) + P2(:))'/2;
    mid2 = (P1(:) + P3(:))'/2;
    mid3 = (P2(:) + P3(:))'/2;
    a_x=[a_x mid1(1) mid2(1) mid3(1)];
    a_y=[a_y mid1(2) mid2(2) mid3(2)];
end

% n_ele=size(a_x);
% for i=1:1:n_ele-1
%    for j=i+1:1:ele
%        if a_x(i)==a_x(j) && a_y(i)==a_y(j)
%    end
% end

a_x2=[a_x2 a_x];
a_y2=[a_y2 a_y];

DT=delaunayTriangulation(a_x2(:),a_y2(:),C);
%triplot(DT)
figure(10)
TF = isInterior(DT);
% 
Tria2 = triangulation(DT.ConnectivityList(TF,:),DT.Points(:,1),DT.Points(:,2));
triplot(Tria2)

end

