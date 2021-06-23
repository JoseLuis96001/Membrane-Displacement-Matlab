function [A,alph_array] = SolveSys(Tria,Point,Adja,Neig,Nmat,min,max,tou)
if tou==1
alpha=0.5;
r=1;
a=-2.3;
b=1.6;
%%% Linear System %%%%%%%%%%%%%%
A=zeros(Nmat,Nmat); %Matrix A
fs=zeros(Nmat,1);%vector f

for i=1:length(Tria.ConnectivityList)
   ik1=Tria(i,1);%index point1
   ik2=Tria(i,2);%index point2
   ik3=Tria(i,3);%index point3
   
   %x y coordinates for the 3 points ik1,...
   v1=Point(ik1,[1 2]);
   v2=Point(ik2,[1 2]);
   v3=Point(ik3,[1 2]);
   
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
   a1=Lkinv(x,y);x1=a1(1);y1=a1(2);
   a2=Lkinv(x,y);x2=a2(1);y2=a2(2);
   a3=Lkinv(x,y);x3=a3(1);y3=a3(2);
   
   p1=p1h(x1,y1);
   p2=p2h(x2,y2);
   p3=p3h(x3,y3);
   
   %gradients
   grad1(xh,yh)=gradient(p1,[x y]);
   grad2(xh,yh)=gradient(p2,[x y]);
   grad3(xh,yh)=gradient(p3,[x y]);
   multi11=dot(grad1(xh,yh),grad1(xh,yh));
   multi22=dot(grad2(xh,yh),grad2(xh,yh));
   multi33=dot(grad3(xh,yh),grad3(xh,yh));
   multi12=dot(grad1(xh,yh),grad2(xh,yh));
   multi13=dot(grad1(xh,yh),grad3(xh,yh));
   multi23=dot(grad2(xh,yh),grad3(xh,yh));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   integral11=int(int(multi11,yh,0,1-xh),xh,0,1);
   integral22=int(int(multi22,yh,0,1-xh),xh,0,1);
   integral33=int(int(multi33,yh,0,1-xh),xh,0,1);
   integral12=int(int(multi12,yh,0,1-xh),xh,0,1);
   integral13=int(int(multi13,yh,0,1-xh),xh,0,1);
   integral23=int(int(multi23,yh,0,1-xh),xh,0,1);
   %triangles inside the boundary
  if((min< ik1)&&(ik1<max) && (min<ik2)&&(ik2<max) && (min<ik3)&&(ik3<max))
   A(ik1-min,ik1-min)=A(ik1-min,ik1-min)+2*area(Tv).*integral11;
   A(ik2-min,ik2-min)=A(ik2-min,ik2-min)+2*area(Tv).*integral22;
   A(ik3-min,ik3-min)=A(ik3-min,ik3-min)+2*area(Tv).*integral33;
   A(ik1-min,ik2-min)=A(ik1-min,ik2-min)+2*area(Tv).*integral12;
   A(ik1-min,ik3-min)=A(ik1-min,ik3-min)+2*area(Tv).*integral13;
   A(ik2-min,ik3-min)=A(ik2-min,ik3-min)+2*area(Tv).*integral23;
   A(ik2-min,ik1-min)=A(ik2-min,ik1-min)+2*area(Tv).*integral12;
   A(ik3-min,ik1-min)=A(ik3-min,ik1-min)+2*area(Tv).*integral13;
   A(ik3-min,ik2-min)=A(ik3-min,ik2-min)+2*area(Tv).*integral23;
  end
  %triangles with 2 points in boundary
  for j=1:3
      if(Neig(i,j)== -1 && j==1)
          A(ik3-min,ik3-min)=A(ik3-min,ik3-min)+2*area(Tv).*integral33;
      end  
      if(Neig(i,j)== -1 && j==2)
          A(ik1-min,ik1-min)=A(ik1-min,ik1-min)+2*area(Tv).*integral11;
      end
      if(Neig(i,j)== -1 && j==3)
          A(ik2-min,ik2-min)=A(ik2-min,ik2-min)+2*area(Tv).*integral22;
      end
  end
  %triangles with one point in the boundary
  if((ik1<=min || ik1>=max) &&(min<ik2)&&(ik2<max) && (min<ik3)&&(ik3<max) )
   A(ik2-min,ik2-min)=A(ik2-min,ik2-min)+2*area(Tv).*integral22;
   A(ik3-min,ik3-min)=A(ik3-min,ik3-min)+2*area(Tv).*integral33;
   A(ik2-min,ik3-min)=A(ik2-min,ik3-min)+2*area(Tv).*integral23;
   A(ik3-min,ik2-min)=A(ik3-min,ik2-min)+2*area(Tv).*integral23;
  end 
   if((ik2<=min || ik2>=max) &&(min< ik1)&&(ik1<max)&& (min<ik3)&&(ik3<max) )
   A(ik1-min,ik1-min)=A(ik1-min,ik1-min)+2*area(Tv).*integral11;
   A(ik3-min,ik3-min)=A(ik3-min,ik3-min)+2*area(Tv).*integral33;
   A(ik1-min,ik3-min)=A(ik1-min,ik3-min)+2*area(Tv).*integral13;
   A(ik3-min,ik1-min)=A(ik3-min,ik1-min)+2*area(Tv).*integral13;
   end 
   if((ik3<=min || ik3>=max) && (min< ik1)&&(ik1<max)&&(min<ik2)&&(ik2<max))
   A(ik1-min,ik1-min)=A(ik1-min,ik1-min)+2*area(Tv).*integral11;
   A(ik2-min,ik2-min)=A(ik2-min,ik2-min)+2*area(Tv).*integral22;
   A(ik1-min,ik2-min)=A(ik1-min,ik2-min)+2*area(Tv).*integral12;
   A(ik2-min,ik1-min)=A(ik2-min,ik1-min)+2*area(Tv).*integral12;
 
   end 
  
   
   %Funcion de fuerza 
for j=1:3
    %global numerated nodes inside the boundary 
    if (min< G_i(j,1))&&(G_i(j,1)<max)
    if ((G(j,1)-a)^2+(G(j,2)-b)^2 <= r^2)
        p1s(xh,yh) = 1 - xh - yh; %Pol.1.sombrero
        p2s(xh,yh) = xh; %Pol.2.sombrero
        p3s(xh,yh) = yh; %Pol.3.sombrero
        
        f_args = Lk(xh,yh);
        f_xh = f_args(1);
        f_yh = f_args(2);
        f(x,y) = alpha*(cos( (pi/2)*(((x-a)^2+(y-b)^2)/(r^2))).^2);
        inte1 = f(f_xh,f_yh);
        intepol1 = p1s(xh,yh);
        intepol2 = p2s(xh,yh);
        intepol3 = p3s(xh,yh);
        integrandofpol1(xh,yh) = f(f_xh,f_yh)*p1s(xh,yh);
        integrandofpol2(xh,yh) = f(f_xh,f_yh)*p2s(xh,yh);
        integrandofpol3(xh,yh) = f(f_xh,f_yh)*p3s(xh,yh);
        %Definir valores de referencia
        Valores_referencia = [0,0;1,0;0,1;0,0.5;0.5,0.5;0.5,0;0.5,0;1/3,1/3]; %Valores de la componente
        Escalar = [3;3;3;8;8;8;27]; %Escalares de la suma
        IntegralF1=0;
        IntegralF2=0;
        IntegralF3=0;
        for e=1:7
            point = Valores_referencia(e,:);
            valor1 = integrandofpol1(point(1,1),point(1,2))*Escalar(e,1)*(1/120);
            IntegralF1 = IntegralF1 + valor1;
            valor2 = integrandofpol2(point(1,1),point(1,2))*Escalar(e,1)*(1/120);
            IntegralF2 = IntegralF2 + valor2;
            point = Valores_referencia(e,:);
            valor3 = integrandofpol3(point(1,1),point(1,2))*Escalar(e,1)*(1/120);
            IntegralF3 = IntegralF3 + valor3; 
        end
        IntegralFinal = IntegralF1 + IntegralF2 + IntegralF3; 
    else
        IntegralFinal = 0; 
    end
    CalculatePoint = G_i(j,1);
    fs(CalculatePoint-min) = fs(CalculatePoint-min) + 2.*area(Tv).*IntegralFinal; 
    
    end  
end
   
   
   
end

%%%% Linear system%%%%%%%%%%%%%%%%%
% Ax = B
alpha_array = A\fs; %Gauss



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alph_array=alpha_array';
%alph_array with the complete dimension
alph_array=[zeros(1,min) alph_array zeros(1,length(Point)-max+1)];
alph_array=alph_array'
% trisurf(Tria.ConnectivityList(), Point(:,1), Point(:,2),alph_array')
% colorbar;  
% colormap jet;
end
end