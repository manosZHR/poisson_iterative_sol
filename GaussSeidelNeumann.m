function [u,res,iter]= GaussSeidelNew(Nr,Nt,rorigin,rlast,torigin,tlast,itgmr,omega,tol)

    ii=[]; ismax=[]; x=[]; y =[]; u=[]; rhs=[]; qrr=[]; 
    du=[];


    string();
    
    %GRID Coordinates

    dr =  (rlast-rorigin)/(Nr-1);
    dt = (tlast-torigin)/(Nt-1);

    for i=2:Nt+1
        for j=2:Nr+1
            L=ii(i,j);
            T(L)=torigin+dt*(i-2);
            R(L)=rorigin+dr*(j-2);
            u(L)=0;
        end
    end
   
    
    rhsDiri();

    solver();

    %% String Numbering

        function string()
    
            %initialization of ii
            for i=1:Nt+2
                for j=1:Nr+2
                    ii(i,j) = 1;
                end
            end
    
            ismax=0; %1D counter of nodes
            for i=2:Nt+1
                for j=2:Nr+1
                    ismax=ismax+1;
                    ii(i,j)=ismax;
                end
            end
        end

%% Boundary Conditions
   
        function rhsDiri()

            for i=2:Nt+1
                for j=2:Nr+1 
                    L=ii(i,j);
                        if  j==2  %Dirichlet
                            rhs(L)=R(L)^2;
                            u(L)=rhs(L);

                        elseif j==Nr+1 %Neumann
                            rhs(L)=0;

                        else
                            rhs(L)=-4;
                        end
                end
            end
    
          
        end

%% Solver
%A'Du = b - Au
%Calculation of norm of b-Au (residual)
%Calculation of u

function solver()

    itmmn = 1;

    for L=1:ismax
        qrr(L) = 0;  %initialization
    end

    for iter=itmmn:itgmr

            %Delta Formulation - Compute Residual  %Au-b
            resit=0;

           for i=2:Nt+1
                for j=2:Nr+1

                    L=ii(i,j);
                    ie=ii(i+1,j);
                    iw=ii(i-1,j);
                    in=ii(i,j+1);
                    is=ii(i,j-1);

                    ec=-2*( 1/dr^2 + 1/((R(L)^2)*(dt^2)) );
                    hc=1/((R(L)^2)*(dt^2));
                    bc=1/((R(L)^2)*(dt^2));
                    fc=1/dr^2+1/(R(L)*2*dr);
                    dc=1/dr^2-1/(R(L)*2*dr);
                    
                    if  j==2 
                        L=ii(i,j);
                        
                        qrr(L) = u(L)-rhs(L);

                    elseif j==Nr+1
                        
                        L=ii(i,j);
                        is=ii(i,j-1);
                        iss=ii(i,j-2);
        
                        c_i = 3/(dr*2);
                        c_is = -2/dr;
                        c_iss = 1/(2*dr);

                        qrr(L) = -( ...
                            ...
                            c_i*u(L)+c_is*u(is)+c_iss*u(iss)+rhs(L));

                    elseif i==Nt+1 & j ~= 2 & j ~= Nr+1 

                        L=ii(i,j);
                        ie=ii(3,j);
                        iw=ii(i-1,j);
                        in=ii(i,j+1);
                        is=ii(i,j-1);

                        qrr(L) = -( ...
                            ...
                            hc*u(ie)+bc*u(iw)+fc*u(in)+dc*u(is)  ... 
                            +ec*u(L)+rhs(L));

                    elseif i==2 & j ~= 2 & j ~= Nr+1 

                        L=ii(i,j);
                        ie=ii(i+1,j);
                        iw=ii(Nt,j);
                        in=ii(i,j+1);
                        is=ii(i,j-1);

                        qrr(L) = -( ...
                            ...
                            hc*u(ie)+bc*u(iw)+fc*u(in)+dc*u(is)  ... 
                            +ec*u(L)+rhs(L));
                    else

                        qrr(L) = -( ...
                            ...
                            hc*u(ie)+bc*u(iw)+fc*u(in)+dc*u(is)  ... 
                            +ec*u(L)+rhs(L));
                    end

                    resit=resit+qrr(L)*qrr(L);
                    
                end

           end

           resit=sqrt(resit)/ismax;

           res(iter)=resit;

           if resit<tol; break; end

            GaussSeidel()

            for L=1:ismax
                u(L)=u(L)+du(L);
            end


    end


end

    function GaussSeidel()

        for L=1:ismax+1
            du(L)= 0; %du initialization
        end

        %A' calculation for delta formulation 
         

            for i=2:Nt+1
                for j=2:Nr+1   
                    L=ii(i,j);
                    ie=ii(i+1,j);
                    iw=ii(i-1,j);
                    in=ii(i,j+1);
                    is=ii(i,j-1);
                   
                    %capital letters
                    if  j==2    %dirichlet boundary
                        ec=1;
                        hc=0;
                        bc=0;
                        fc=0;
                        dc=0;

                    elseif j==Nr+1   %neumann boundary
                        ec = 1/dr;
                        dc = -1/dr;
                        hc=0;
                        bc=0;
                        fc=0;


                    elseif i==Nt+1 & j ~= 2 & j ~= Nr+1  %split line - last node of the circle


                        ec=-2*( 1/dr^2 + 1/((R(L)^2)*(dt^2)) );
                        hc=0;
                        bc=1/((R(L)^2)*(dt^2));
                        fc=1/dr^2+1/(R(L)*2*dr);
                        dc=1/dr^2-1/(R(L)*2*dr);

                    elseif i==2 & j ~= 2 & j ~= Nr+1   %split line - first node of the circle

                        ec=-2*( 1/dr^2 + 1/((R(L)^2)*(dt^2)) );
                        hc=1/((R(L)^2)*(dt^2));
                        bc=0;
                        fc=1/dr^2+1/(R(L)*2*dr);
                        dc=1/dr^2-1/(R(L)*2*dr);

                    else  %internal nodes
                        ec=-2*( 1/dr^2 + 1/((R(L)^2)*(dt^2)) );
                        hc=1/((R(L)^2)*(dt^2));
                        bc=1/((R(L)^2)*(dt^2));
                        fc=1/dr^2+1/(R(L)*2*dr);
                        dc=1/dr^2-1/(R(L)*2*dr);
                    end

                    du(L) = (1/ec)*(qrr(L)-hc*du(ie)-bc*du(iw)-fc*du(in)-dc*du(is));

                end

            end



    end


end