function [Uhd] = evaluation(Punto,Tria,Neig,u,h)
%Valores
    Xo=Punto(1,1);
    Yo=Punto(1,2);
    %ID=pointLocation(Tria,Xo,Yo);
    ID=findpoint(Tria,Neig,Punto,h);
    ik1=Tria(ID,1);%index point1
    ik2=Tria(ID,2);%index point2
    ik3=Tria(ID,3);%index point3
    
   %x y coordinates for the 3 points ik1,...
   v1=Tria.Points(ik1,[1 2]);
   v2=Tria.Points(ik2,[1 2]);
   v3=Tria.Points(ik3,[1 2]);
   
   G=[v1;v2;v3]; %matrix with points k1,k2,k3 
   G_i=[ik1;ik2;ik3]; %matrix of global numerated nodes
   
   Tv=polyshape(G(:,1),G(:,2));
   AreaTi=area(Tv);%Area of Tria(k)
   
   %Base on P1(T^)
   syms xh yh x y p1 p2 f
   TR = [0 0;0 1;1 1];%Reference triangle
   p1h(xh,yh)=1-xh-yh; %p1^
   p2h(xh,yh)=xh; %p2^
   p3h(xh,yh)=yh; %p3^
   
   %Linear application Lk
   Ak=[v2(1,1)-v1(1,1), v3(1,1)-v1(1,1); v2(1,2)-v1(1,2), v3(1,2)-v1(1,2)];
   Xh=[xh;yh];
   Vb=[v1(1,1);v1(1,2)];
   Lk(xh,yh)=Ak*Xh+Vb;
   %Inverse of Lk
   %Delta=1/det(Ak);
   X=[x-v1(1,1);y-v1(1,2)];
   Lkinv(x,y)=inv(Ak)*X;
   %Base on T satisfying the interpolation conditions
   a1=Lkinv(Xo,Yo);
   a2=Lkinv(Xo,Yo);
   a3=Lkinv(Xo,Yo);
   

%    
   x1=a1(1);y1=a1(2);
   x2=a2(1);y2=a2(2);
   x3=a3(1);y3=a3(2);
   
   p1=p1h(x1,y1);
   p2=p2h(x2,y2);
   p3=p3h(x3,y3);
   
  
    
    Uh = u(ik1,1)*p1 + u(ik2,1)*p2 + u(ik3,1)*p3; 
    Uhd=double(Uh);
    
end
    
    
    
    
    
    
    