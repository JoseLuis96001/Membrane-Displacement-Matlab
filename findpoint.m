function [To] = findpoint(Tria,Neig,Point,h)
%clf(10);
% [f,~]=size(Tria);
% To=randi([1 f]);
%To=1;%74;
l=h;
[f,~]=size(Tria);
seg=[];
for i=1:1:f
    vals=Tria(i,:);
    c=0;
    for j=1:1:3
        if vals(j)==1 || vals(j)==2
            c=c+1;
        elseif vals(j)==1+l || vals(j)==2+l
            c=c+1;
        elseif vals(j)==1+l*2 || vals(j)==2+l*2
            c=c+1;
        end
    end
    if c==2
        seg=[seg i];
    end
end

if (Point(1)<0 && Point(2)>=0) || (Point(1)<0 && Point(2)<0) 
    To=seg(2);
elseif Point(1)>=0 && Point(2)>=0
    To=seg(3);
elseif Point(1)>=0 && Point(2)<0
    To=seg(1);
end
Li=[];
while 1
    p1=Tria(To,1);
    p2=Tria(To,2);
    p3=Tria(To,3);
    xy_p1=Tria.Points(p1,:);
    xy_p2=Tria.Points(p2,:);
    xy_p3=Tria.Points(p3,:);
    figure(10);
    triplot(Tria);
    hold on
    loc=incenter(Tria,To);
    plot(loc(1),loc(2),'r*');
    display(To);
    Li=[Li To];
    pause(1)
    reep1 = isPositive(Point,xy_p2,xy_p3);
    reep2 = isPositive(xy_p1,Point,xy_p3);
    reep3 = isPositive(xy_p1,xy_p2,Point);

    if reep1==true && reep2==true && reep3==true
        break
    end

    if reep1==false  
        if Neig(To,2) ~= -1
%             To=Neig(To,2)
            p=Neig(To,2);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end    
            To=p;
            display(["reep1 :"+num2str(To)])
        elseif  Neig(To,1) ~= -1
            p=Neig(To,1);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end    
            To=p
        elseif  Neig(To,3) ~= -1
            p=Neig(To,3);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end
            To=p
        end
    elseif reep2==false  
        if Neig(To,3) ~= -1
            %To=Neig(To,3)
            p=Neig(To,3);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end    
            To=p;
            display(["reep2 :"+num2str(To)])
        elseif  Neig(To,2) ~= -1
            p=Neig(To,2);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end    
            To=p
        elseif  Neig(To,1) ~= -1
            p=Neig(To,1);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end
            To=p
        end
    elseif reep3==false 
        if Neig(To,1) ~= -1
%             To=Neig(To,1)
            p=Neig(To,1);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end    
            To=p;
            display(["reep3 :"+num2str(To)])
        elseif  Neig(To,2) ~= -1
            p=Neig(To,2);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end    
            To=p
        elseif  Neig(To,3) ~= -1
            p=Neig(To,3);
            if Neig(p,1)==To
                Neig(p,1)=-1;
            elseif Neig(p,2)==To
                Neig(p,2)=-1;
            elseif Neig(p,3)==To
                Neig(p,3)=-1;
            end
            To=p
        end    
    end
end 
%end fun
end