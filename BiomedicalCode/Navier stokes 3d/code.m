clear
clear
M = 10;
N =10;  
D=0.04; L=40;
eta=1.85e-5;
u0=5;
P1=20;
P2=1;
starte = 12;
longeur= 50;
dx=L/(N-1); dy=D/(M-1);
K=21;
F=zeros(4*M*K*N,1);
[A,B]=matrice3dd(F,M,N,K,L,D,P1,P2);

%[A,B]=deformerr(A,B,M,N,L,D,P1,P2,starte,pas+2,longeur);  
%[A,B]=deformeeerr(A,B,M,N,L,D,P1,P2,3,3,8); 
%[A,B]=deformerr(A,B,M,N,L,D,P1,P2,34,3,10); 
%[A,B]=deformerr(A,B,M,N,L,D,P1,P2,50,3,10); 

%[A,B]=deformeer(A,B,M,N,L,D,P1,P2,starte+2,pas,longeur-5);  
%[A,B]=deformeer(A,B,M,N,L,D,P1,P2,24,5,20);  
%[A,B]=deformeer(A,B,M,N,L,D,P1,P2,40,3,8); 
%[A,B]=deformeer(A,B,M,N,L,D,P1,P2,2,2,6);  
%[A,B]=deformeer(A,B,M,N,L,D,P1,P2,8,4,14);  
t=0;
for j=1:4*M*K*N
    if A(:,j)==0
        disp(j);
    end
    if A(j,:)==0
        disp(j);
    end
end
S=A\B;
Ur = zeros(M*K*N,1);
for j=3*M*K*N+1:4*M*K*N
    Ur(j-3*M*K*N)=S(j);
end
pi=3.14;
% x=zeros(N*M*K,1);
% y=zeros(N*M*K,1);
% z=zeros(N*M*K,1);
X=[];Y=[];Z=[];
pp=linspace(1,L,M);
pi=3.14;
for k=1:M
    for j =1:K
        for i=1:N
             r=(i-1)*D/2*N;
             theta=(j-1)*2*pi/(K-1);
             x = r*cos(theta);
             y = r*sin(theta);
             z = pp(k);
             X = [X, x];
             Y = [Y, y];
             Z = [Z, z];
        end
    end
end
scatter3(X,Y,Z,10,Ur)
%U=S(1:M*N,1);
%U1=reshape(U,M,N);
%P=S(2*M*N+1:3*M*N,1);
%PP=reshape(P,M,N);
%V=S(M*N+1:2*M*N,1);
%V1=reshape(V,M,N);
%C=sqrt(U.^2 + V.^2);
%C1=reshape(C,M,N);
%surf(X,Y,U1); title('Vitesse ux');  colorbar;
%[X,Y] = meshgrid(1:	M,1:N);

%quiver(X,Y,U1,V1,0)
%pcolor(X,Y,U1);