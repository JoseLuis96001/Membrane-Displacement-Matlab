%%%%%%%%%%%%%%%%%Complete process and correct order of execution %%%%%%%%%%%%%%%%%%%
%%Initial data
R=5;
l=2.5;
h=5;
%%GenMesh
%GenMesh returns:
%Structures: Tria Point Adja Neig
%integer values Nmat min max 
%Constrain matrix C
%value tou true if all elements in Tria are positive oriented
[Tria,Point,Adja,Neig,Nmat,min,max,C,tou]=GenMesh(R,l,h);%1,1,4 2,1,5,

%%Solve the system
%It returns the matrix A and the alphas found
[A,alph_array] = SolveSys(Tria,Point,Adja,Neig,Nmat,min,max,tou);

%%some Graphs
figure(3)
spy(A) %Matrix A 
figure(4)
trisurf(Tria.ConnectivityList(), Point(:,1), Point(:,2),alph_array) %3D Graph
xlabel('x'), ylabel('y'), zlabel('z')
colorbar;  
colormap jet;

  
%%Evaluation
a=-2.3;
b=1.6;
Ptoval=[a b];%this is the point to evaluate
[Uh] = evaluation(Ptoval,Tria,Neig,alph_array,h);


