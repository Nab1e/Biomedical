function [A,B]=matrice2d(F,nx,ny,L,D,P1,P2)

dx=L/(ny-1);
dy=D/(nx-1);
%disp(dx);
%disp(dy);
eta=0.0000185;
N=3*nx*ny;
A=spalloc(N,N,7*N);
%A=zeros(N,N);
B=reshape(F,N,1);
%disp(eta/dx^2);
g=9.8;

coeff=[eta/dx^2,eta/dy^2,(-2/dx^2-2/dy^2)*eta,eta/dy^2,eta/dx^2,1/(2*dx),-1/(2*dx)];
num = [  -nx, -1, 0, 1, nx, -nx+2*nx*ny, nx+2*nx*ny ];

coeff2=[eta/dx^2,eta/dy^2,(-2/dx^2-2/dy^2)*eta,eta/dy^2,eta/dx^2,1/(2*dy),-1/(2*dy)];
num2 = [  -nx, -1, 0, 1, nx, -1+nx*ny, 1+nx*ny ];

    for j=2:ny-1
        k=(j-1)*nx+2;
        %disp(k);
       % A((j-1)*nx+nx+2*nx*ny,k+num)=coeff;
        k2=k+nx*ny;
        %disp(k2);
        A((j-1)*nx+1+2*nx*ny,k2+num2)=coeff2;
   end


% C.L. sur les frontieres
for i=2:nx-1
    k=i; A(k,:)=0; A(k,k)=-1.0; A(k,k+nx)=1.0; A(k,(ny-1)*nx+i)=-1.0; A(k,(ny-2)*nx+i)=1.0; B(k)=0; %derivee de ux nulle
    k=i+nx*ny; A(k,:)=0;  A(k,k)=-1.0; A(k,k+nx)=1.0; A(k,(ny-1)*nx+k)=-1.0; A(k,(ny-2)*nx+k)=1.0; B(k)=0;%derivee de uy nulle
    k=i+2*nx*ny;  A(k,i)=1.0; A(k,(ny-1)*nx+i)=-1.0; B(k)=0;   %ux=ux initiale
    k=(ny-1)*nx+i;   A(k,i+nx*ny)=1.0;A(k,(ny-1)*nx+i+nx*ny)=-1.0; B(k)=0; % uy=uy
    k=(ny-1)*nx+i+nx*ny;   A(k,i+2*nx*ny)=1.0; B(k)= P1;  %P=P1
    k=(ny-1)*nx+i+2*nx*ny;  A(k,k)=1.0; B(k)= P2;    %P=P2
end
    k=(ny-1)*nx+1+2*nx*ny; A(k,k)=1.0; B(k)= P2;    %P=P2
    k=(ny-1)*nx+nx+2*nx*ny;   A(k,k)=1.0; B(k)= P2;  %P=P2
    k=1+2*nx*ny;  A(k,k)=1.0; B(k)= P1;    %P=P1
    k=nx+2*nx*ny;   A(k,k)=1.0; B(k)= P1;  %P=P1

% conditions aux limites ux=0 uy=0 dP/dx=cte en haut et en bas du tube


for j=1:ny
    k=(j-1)*nx+1; A(k,:)=0; A(k,k)=1.0; B(k)=0;%ux
    k=(j-1)*nx+1+nx*ny; A(k,:)=0; A(k,k)=1.0; B(k)=0;%uy
    k=(j-1)*nx+nx; A(k,:)=0; A(k,k)=1.0; B(k)=0;%ux
    k=(j-1)*nx+nx+nx*ny; A(k,:)=0; A(k,k)=1.0; B(k)=0;%uy
   
end


 for j=2:ny-1
      %k=(j-1)*nx+1+2*nx*ny; A(k,:)=0; A(k,(j-1)*nx+2+nx*ny)=1.0; A(k,(j-1)*nx+1+nx*ny)=-1.0; B(k)=0;
      k=(j-1)*nx+nx+2*nx*ny; A(k,:)=0; A(k,(j-1)*nx+nx+nx*ny)=1.0; A(k,(j-1)*nx+nx-1+nx*ny)=-1.0; B(k)=0;
 end
 %j=1; k=(j-1)*nx+1; A(k,:)=0; A(k,(j-1)*nx+2+nx*ny)=1.0; A(k,(j-1)*nx+1+nx*ny)=-1.0; B(k)=0;
% j=ny; k=(j-1)*nx+1; A(k,:)=0; A(k,(j-1)*nx+2+nx*ny)=1.0; A(k,(j-1)*nx+1+nx*ny)=-1.0; B(k)=0;

      k=nx+nx-1;
A(nx+nx+2*nx*ny,k+num)=coeff;
%    
%
%
coeff=[eta/dx^2,eta/dy^2,(-2/dx^2-2/dy^2)*eta,eta/dy^2,eta/dx^2,1/(dx),-1/(dx)];
num = [  -nx, -1, 0, 1, nx,  2*nx*ny, nx+2*nx*ny ];

coeff2=[eta/dx^2,eta/dy^2,(-2/dx^2-2/dy^2)*eta,eta/dy^2,eta/dx^2,1/(dy),-1/(dy)];
num2 = [  -nx, -1, 0, 1, nx,  nx*ny, 1+nx*ny ];

coeff3=[-1/(dx), 1/(dx), -1/(dy), 1/(dy)];
num3 = [  -nx-2*nx*ny, nx-2*nx*ny, -1-nx*ny, 1-nx*ny ];




for i=2:nx-1
    for j=2:ny-1
        k=(j-1)*nx+i;
        %disp(k);
        A(k,k+num)=coeff;
        B(k)=0;
        k2=k+nx*ny;
        %disp(k2);
        A(k2,k2+num2)=coeff2;
        B(k2)=0;
        k3=k2+nx*ny;
        %disp(k3);
        A(k3,k3+num3)=coeff3;  
    end
end