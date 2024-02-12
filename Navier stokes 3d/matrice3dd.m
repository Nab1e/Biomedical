function [A,B]=matrice3dd(F,nk,nr,nt,L,D,P1,P2)


PI=3.14;


y = linspace(1,0,nk);


dr=D/(nr-1);


dk=L/(nk-1);


dtheta=PI/(nt-1);


 

 

 

N=4*nr*nt*nk;


eta=0.0000185;


A=spalloc(N,N,7*N);

B=reshape(F,N,1);


% C.L. sur les frontieres





 

 

 

 

for i=2:nr-1


    for j=1:nt


        for k=2:nk-1

           

            r=(i-1)*D/2*nr;

           
            theta=(j-1)*2*pi/(nt-1); 

           

            coeff1=[1/dk^2,1/(r^2)*dtheta^2,1/(2*r*dr)+1/(dr^2),-2*(1/(dk^2)+1/(dr^2)+1/(dtheta^2)),-1/(2*r*dr)+1/(dr^2),-1/((r^2)*dtheta^2),1/(dk^2),-1/2*eta*dr,1/2*eta*dr];


            num1 = [nr*nt,nr,1,0,-1,-nr,-nr*nt,1+3*nr*nt*nk,-1+3*nr*nt*nk];


            Phi=(k-1)*nt*nr + (j-1)*nr + i;


            A(Phi,Phi+num1)=coeff1;


            

            coeff2 = [1/dk^2,1/(r^2)*dtheta^2,1/(2*r*dr)+1/(dr^2),-2*(1/(dk^2)+1/(dr^2)+1/(dtheta^2)),-1/(2*r*dr)+1/(dr^2),-1/((r^2)*dtheta^2),1/(dk^2),1/(2*r*dtheta),-1/(2*r*dtheta)];


            num2 = [nr*nt + nr*nt*nk, nr + nr*nt*nk, 1+nr*nt*nk, 0, -1+nr*nt*nk,-nr+nr*nt*nk, -nr*nt+nr*nt*nk,3*nr*nt*nk-nr,3*nr*nt*nk+nr];


       

            Phi2=(k-1)*nt*nr + (j-1)*nr + i + nr*nt*nk;


            A(Phi2,Phi+num2)=coeff2;

            

            coeff3=[1/dk^2,1/(r^2)*dtheta^2,1/(2*r*dr)+1/(dr^2),-2*(1/(dk^2)+1/(dr^2)+1/(dtheta^2)),-1/(2*r*dr)+1/(dr^2),-1/((r^2)*dtheta^2),1/(dk^2),-1/(2*eta*dk),1/(2*eta*dk)];


            num3 = [nr*nt+2*nr*nt*nk,nr+2*nr*nt*nk,1+2*nr*nt*nk,2*nr*nt*nk,-1+2*nr*nt*nk,-nr+2*nr*nt*nk,-nr*nt+2*nr*nt*nk,nr*nt+3*nr*nt*nk,-nr*nt+3*nr*nt*nk];


            Phi3= Phi + 2*nr*nt*nk;

            A(Phi3,Phi+num3)=coeff3;

 

 

        end

    end

end


 

 for k=1:nk

    for j=1:nt

        Phi=(k-1)*nt*nr + (j-1)*nr + nr+3*nr*nt*nk; A(Phi,:)=0; A(Phi,Phi)=y(k); B(Phi)=0;%pression


    end

end

for j=1:nt

     for i=1:nr


         Phi=(j-1)*nr +i +3*nr*nt*nk; A(Phi,:)=0; A(Phi,Phi)=1.0 ; B(Phi)=1;


         Phi=(nk-1)*nt*nr + (j-1)*nr +i +3*nr*nt*nk; A(Phi,:)=0;A(Phi,Phi)=1.0 ; B(Phi)=0;
         

     end

end

for j=1:nt

    for k=2:nk-1

         for i=1:nr

          r=(i-1)*D/2*nr;


          theta=(j-1)*2*pi/(nt-1);


         % coeff4=[1/r,1/(2*dr),-1/(2*dr),1/(2*dk),-1/(2*dk)];


          %num4 = [0,1,-1,nr*nt+2*nr*nt*nk,-nr*nt+2*nr*nt*nk];

          coeff4=[1/dr^2,1/dtheta^2,1/dk^2,-2*(1/dr^2 + 1/dtheta^2 + 1/dk^2),1/dk^2,1/dtheta^2, 1/dr^2];

          num4 = [1,nr,nr*nt,0,-nr*nt,-nr, -1];

        Phi=(k-1)*nt*nr + (j-1)*nr + i;

       

        A(Phi+3*nr*nt*nk,Phi+3*nr*nt*nk+num4)=coeff4;

        

       end

        

    end


end


 

for k=1:nk

    for j=1:nt

        i=nr;

        Phi=(k-1)*nt*nr + (j-1)*nr + i;A(Phi,:)=0; A(Phi,Phi)=1.0; B(Phi)=0;%ur

        Phi=(k-1)*nt*nr + (j-1)*nr + i + nr*nt*nk;A(Phi,:)=0; A(Phi,Phi)=1.0; B(Phi)=0;%ut

        Phi=(k-1)*nt*nr + (j-1)*nr + i + 2*nr*nt*nk;A(Phi,:)=0; A(Phi,Phi)=1.0; B(Phi)=0;%uz
        
        i=1;
        
        Phi=(k-1)*nt*nr + (j-1)*nr + i + nr*nt*nk;A(Phi,:)=0; A(Phi,Phi)=1.0; B(Phi)=0;%ut


    end

end

for k=1:nk

    for j=1:(nt-1)/2

        i=1;

        r=(i-1)*D/2*nr;

        Phi=(k-1)*nt*nr + (j-1)*nr + i; A(Phi,:)=0; A(Phi,Phi)=1.0; A(Phi,Phi + 1 + nr*nt/2-nr/2)=-1.0; B(Phi)=0;%ur

        Phi=(k-1)*nt*nr + (j-1)*nr + i + nt*nr*nt; A(Phi,:)=0; A(Phi,Phi)=1.0; A(Phi,Phi + 1 + nr*nt/2-nr/2)=-1.0; B(Phi)=0;%ur

        Phi=(k-1)*nt*nr + (j-1)*nr + i + 2*nt*nr*nt; A(Phi,:)=0; A(Phi,Phi)=1.0; A(Phi,Phi + 1 + nr*nt/2-nr/2)=-1.0; B(Phi)=0;%ur

        Phi=(k-1)*nt*nr + (j-1)*nr + i + 3*nt*nr*nt; A(Phi,:)=0; A(Phi,Phi)=1.0; A(Phi,Phi + 1 + nr*nt/2-nr/2)=-1.0; B(Phi)=0;%ur


    end

    for j=(nt-1)/2+1:nt

        i=1;

        r=(i-1)*D/2*nr;

        Phi=(k-1)*nt*nr + (j-1)*nr + i; A(Phi,:)=0; A(Phi,Phi)=1.0; A(Phi,Phi + 1 -nr*nt/2+nr/2)=-1.0; B(Phi)=0;%ur

        Phi=(k-1)*nt*nr + (j-1)*nr + i + nt*nr*nt; A(Phi,:)=0; A(Phi,Phi)=1.0; A(Phi,Phi + 1  -nr*nt/2+nr/2)=-1.0; B(Phi)=0;%ur

        Phi=(k-1)*nt*nr + (j-1)*nr + i + 2*nt*nr*nt; A(Phi,:)=0; A(Phi,Phi)=1.0; A(Phi,Phi + 1  -nr*nt/2+nr/2)=-1.0; B(Phi)=0;%ur

        Phi=(k-1)*nt*nr + (j-1)*nr + i + 3*nt*nr*nt; A(Phi,:)=0; A(Phi,Phi)=1.0; A(Phi,Phi + 1  -nr*nt/2+nr/2)=-1.0; B(Phi)=0;%ur


    end

end

for j=1:nt

     for i=1:nr


     %k=1


     %ur


     Phi=(j-1)*nr + i;  A(Phi,:)=0; A(Phi,Phi)=-1.0; A(Phi,(2-1)*nt*nr + (j-1)*nr + i)=1.0; A(Phi,(nk-1)*nt*nr + (j-1)*nr + i)=-1.0; A(Phi,(nk-2)*nt*nr + (j-1)*nr + i)=1.0;B(Phi)=0;


    %ut


     Phi=(j-1)*nr + i +nr*nt*nk;  A(Phi,:)=0; A(Phi,Phi)=-1.0; A(Phi,(2-1)*nt*nr + (j-1)*nr + i+nr*nt*nk)=1.0; A(Phi,(nk-1)*nt*nr + (j-1)*nr + i+nr*nt*nk)=-1.0; A(Phi,(nk-2)*nt*nr + (j-1)*nr + i+nr*nt*nk)=1.0;B(Phi)=0;


    %uz


     Phi=(j-1)*nr + i +2*nr*nt*nk;  A(Phi,:)=0; A(Phi,Phi)=-1.0; A(Phi,(2-1)*nt*nr + (j-1)*nr + i+2*nr*nt*nk)=1.0; A(Phi,(nk-1)*nt*nr + (j-1)*nr + i+2*nr*nt*nk)=-1.0; A(Phi,(nk-2)*nt*nr + (j-1)*nr + i+2*nr*nt*nk)=1.0;B(Phi)=0;


    %pression


 

     %k=Nk


     %ur


     Phi=(nk-1)*nt*nr + (j-1)*nr + i;  A(Phi,:)=0; A(Phi,Phi)=-1.0;A(Phi,(j-1)*nr + i)=1.0;


     %ut


     Phi=(nk-1)*nt*nr + (j-1)*nr + i +nr*nt*nk;  A(Phi,:)=0; A(Phi,Phi)=-1.0;A(Phi,(j-1)*nr + i +nr*nt*nk)=1.0;


     %uz


     Phi=(nk-1)*nt*nr + (j-1)*nr + i  +2*nr*nt*nk;  A(Phi,:)=0; A(Phi,Phi)=-1.0;A(Phi,(j-1)*nr + i +2*nr*nt*nk)=1.0;


      %pression


 

     end

end
end

 