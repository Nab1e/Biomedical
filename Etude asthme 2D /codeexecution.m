clear
M = 60;
N = 60;
DD=zeros(N,1);
RR=zeros(N,1);
DDY=zeros(3*N*M,1);
D=0.018;
L=0.12;
eta=1.85e-5;
u0=5;
P1=1.0822;
P2=1;
Count=0;
deltap=P1-P2;
Resistance=0;
starte = 12;
longeur= 20;
pas=1;
dx=L/(N-1); dy=D/(M-1);
X=[0:dx:L]; Y=[0:dy:D];
theta = -pi/4;
Di=0;
RR=zeros(N,1);
F=zeros(3*M,N);
%
[A,B]=matrice2d(F,M,N,L,D,P1,P2);
%
%[A,B]=deformerr(A,B,M,N,L,D,P1,P2,fix(3*N/4),fix(M/14),fix(N/5));
[A,B]=deformerr(A,B,M,N,L,D,P1,P2,fix(N/8),fix(M/4),fix(N*0.7));
[A,B]=deformerr(A,B,M,N,L,D,P1,P2,fix(N*0.3),fix(M*0.36),fix(N*0.4));
[A,B]=deformerr(A,B,M,N,L,D,P1,P2,fix(N*0.7),fix(M/6),fix(N*0.2));

%[A,B]=deformeer(A,B,M,N,L,D,P1,P2,fix(N*0.4),fix(M*0.36),fix(N*0.6));
[A,B]=deformeer(A,B,M,N,L,D,P1,P2,fix(N*0.5),fix(M/6),fix(2*N/5));
[A,B]=deformeer(A,B,M,N,L,D,P1,P2,fix(N/12),fix(M/4),fix(N*0.7));
%
S=A\B;
U=S(1:M*N,1);
U1=reshape(U,M,N);
P=S(2*M*N+1:3*M*N,1);
PP=reshape(P,M,N);
V=S(M*N+1:2*M*N,1);
V1=reshape(V,M,N);
C=sqrt(U.^2 + V.^2);
C1=reshape(C,M,N);
 surf(X,Y,PP);colorbar
 caxis([P2,P1])
title('Pression'); colorbar;
[X,Y] = meshgrid(1:M,1:N);xlabel(" X ");ylabel(" Y ");
    Count=0;
    for i=1:M
        k=(N-1)*M+i;
        if C(k,1)~=0
           Count=Count+1;
           Vitesse(i,1)=C(k,1);
        end
    end
    Rayonlocale=(Count-1)*dy/2;
    Diametrelocale=2*Rayonlocale;
    debitsurfa=Diametrelocale*sum(Vitesse)/Count;


%     Resistance=Resistance+Resistancelocale;
% end
% K=1:1:N;
% plot(K,RR);
% disp("courbe ");

tot=deltap/mean(debitsurfa);
disp(tot);
disp(mean(debitsurfa));
%  disp("Rapport des resistances "+ tot/resistance);
% % disp("somme des resistances : "+Resistance);
% diff=Resistance-tot;
% %disp("leur différence est " +diff);
% erreur=2*diff/(tot+Resistance);
% disp("erreur est " +erreur);
% disp("************");
% disp(" size de A est "+size(A));
% disp(" size de S est "+size(S));
% DDY=A*S;
% disp(" size de a*s est "+size(DDY));
% disp("resultat de comparaison "+ isequal(DDY,B));
% difference = (DDY-B).^2;
% erreur = mean(difference(:));
% disp("ERREUR entre les deux vecteur colonne "+erreur);
% disp("moyenne de la vitesse dans la trachée est "+mean(C));
% disp("************");