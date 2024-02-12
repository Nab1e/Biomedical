function [A,B]=deformerr(A,B,nx,ny,L,D,P1,P2,starte,pas,longeur)
eta =1.85e-5;

dx=L/(ny-1); dy=D/(nx-1);
coeff=[eta/dx^2,eta/dy^2,(-2/dx^2-2/dy^2)*eta,eta/dy^2,eta/dx^2,1/(2*dx),-1/(2*dx)];
num = [  -nx, -1, 0, 1, nx, -nx+2*nx*ny, nx+2*nx*ny ];

coeff2=[eta/dx^2,eta/dy^2,(-2/dx^2-2/dy^2)*eta,eta/dy^2,eta/dx^2,1/(2*dy),-1/(2*dy)];
num2 = [  -nx, -1, 0, 1, nx, -1+nx*ny, 1+nx*ny ];

qq=0;
fx=1;
fy=longeur;
start=starte;
H=zeros(fx*fy,2);
for v=1:fx
   for w=1:fy
    H(w+qq*fy,1)=v;
    H(w+qq*fy,2)=start+w;
   end  
   qq=qq+1;
end
%disp(H);

for v=1:fx*fy
   i=H(v,1);
   j=H(v,2);
   k=(j-1)*nx+i;
   A(k,:)=0; A(k,k)=1; B(k)=0;
   A(k+nx*ny,:)=0;A(k+nx*ny,k+nx*ny)=1; B(k+nx*ny)=0;
   k=(j-1)*nx+i+2*nx*ny; A(k,:)=0; A(k,k)=1; B(k)=0;
   
   i=H(v,1);
   j=H(v,2)+1;
   k=(j-1)*nx+i;
   A(k,:)=0; A(k,k)=1; B(k)=0;
   A(k+nx*ny,:)=0;A(k+nx*ny,k+nx*ny)=1; B(k+nx*ny)=0;
   
   i=H(v,1);
   j=H(v,2)-1;
   k=(j-1)*nx+i;
   A(k,:)=0; A(k,k)=1; B(k)=0;
   A(k+nx*ny,:)=0;A(k+nx*ny,k+nx*ny)=1; B(k+nx*ny)=0;
   
   i=H(v,1)+1;
   j=H(v,2);
   k=(j-1)*nx+i;
   A(k,:)=0; A(k,k)=1; B(k)=0;
   A(k+nx*ny,:)=0;A(k+nx*ny,k+nx*ny)=1; B(k+nx*ny)=0;
   k=(j-1)*nx+i+2*nx*ny; A(k,:)=0; A(k,(j-1)*nx+i+nx*ny)=1.0; A(k,(j-1)*nx+i-1+nx*ny)=-1.0; B(k)=0;
   
   k=(j-1)*nx+i+1;
   k2=k+nx*ny;
   A((j-1)*nx+i+2*nx*ny,k2+num2)=coeff2;
end

for t=1:pas
qq=0;
fx=t;
fy=longeur-t*2;
start=starte+t;
H=zeros(fx*fy,2);
for v=1:fx
   for w=1:fy
    H(w+qq*fy,1)=v;
    H(w+qq*fy,2)=start+w;
   end  
   qq=qq+1;
end
%disp(H);

for v=1:fx*fy
   i=H(v,1);
   j=H(v,2);
   k=(j-1)*nx+i;
   A(k,:)=0; A(k,k)=1; B(k)=0;
   A(k+nx*ny,:)=0;A(k+nx*ny,k+nx*ny)=1; B(k+nx*ny)=0;
   k=(j-1)*nx+i+2*nx*ny; A(k,:)=0; A(k,k)=1; B(k)=0;
   
   i=H(v,1);
   j=H(v,2)+1;
   k=(j-1)*nx+i;
   A(k,:)=0; A(k,k)=1; B(k)=0;
   A(k+nx*ny,:)=0;A(k+nx*ny,k+nx*ny)=1; B(k+nx*ny)=0;
   
   i=H(v,1);
   j=H(v,2)-1;
   k=(j-1)*nx+i;
   A(k,:)=0; A(k,k)=1; B(k)=0;
   A(k+nx*ny,:)=0;A(k+nx*ny,k+nx*ny)=1; B(k+nx*ny)=0;
   
   i=H(v,1)+1;
   j=H(v,2);
   k=(j-1)*nx+i;
   A(k,:)=0; A(k,k)=1; B(k)=0;
   A(k+nx*ny,:)=0;A(k+nx*ny,k+nx*ny)=1; B(k+nx*ny)=0;
   k=(j-1)*nx+i+2*nx*ny; A(k,:)=0; A(k,(j-1)*nx+i+nx*ny)=1.0; A(k,(j-1)*nx+i-1+nx*ny)=-1.0; B(k)=0;
   
   k=(j-1)*nx+i+1;
   k2=k+nx*ny;
   A((j-1)*nx+i+2*nx*ny,k2+num2)=coeff2;
end
end

end