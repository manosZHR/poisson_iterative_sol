function [u,res,iter]= MSIP9PolarDirichlet(Nr,Nt,rorigin,rlast,torigin,tlast,psi,itgmr,tol)

    ii=[]; ismax=[]; x=[]; y =[]; u=[]; rhs=[]; qrr=[]; 
    am=[]; bm=[]; cm=[]; dm=[]; em=[]; fm=[]; gm=[];  hm=[]; km=[]; aux=[]; du=[];
    
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
    
    

    msip5();
    
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
                        if  j==2 | j==Nr+1  %Dirichlet
                            rhs(L)=R(L)^2;
                            u(L)=rhs(L);

                        else
                            rhs(L)=-4;
                        end
                end
            end
    
          
        end

           
%% SIP Small Letters
%Calculation of SIP Small Letters. Elements of 5-diagonal A' (easily
%invertible version of A)

        function msip5()

         
            for L=1:ismax+1
                am(L)=0;
                bm(L)=0;
                cm(L)=0;
                dm(L)=0;
                em(L)=0;
                fm(L)=0;
                gm(L)=0;
                hm(L)=0;
                km(L)=0;
            end

            %MSIP, stencil 9 nodes
            %A' calculation for delta formulation 
         

            for i=2:Nt+1
                for j=2:Nr+1   
                    L=ii(i,j);
                    ie=ii(i+1,j);
                    iw=ii(i-1,j);
                    in=ii(i,j+1);
                    is=ii(i,j-1);
                    ine=ii(i+1,j+1);
                    ise=ii(i+1,j-1);
                    inw=ii(i-1,j+1);
                    isw=ii(i-1,j-1);
                   
                    %capital letters
                    if  j==2 | j==Nr+1   %dirichlet boundary
                        ec=1;
                        hc=0;
                        bc=0;
                        fc=0;
                        dc=0;
                        kc=0;
                        cc=0;
                        gc=0;
                        ac=0;


                    elseif i==Nt+1 & j ~= 2 & j ~= Nr+1  %split line - last node of the circle


                        ec=-2*( 1/dr^2 + 1/((R(L)^2)*(dt^2)) );
                        hc=0;
                        bc=1/((R(L)^2)*(dt^2));
                        fc=1/dr^2+1/(R(L)*2*dr);
                        dc=1/dr^2-1/(R(L)*2*dr);
                        kc=0;
                        cc=0;
                        gc=0;
                        ac=0;

                    elseif i==2 & j ~= 2 & j ~= Nr+1   %split line - first node of the circle

                        ec=-2*( 1/dr^2 + 1/((R(L)^2)*(dt^2)) );
                        hc=1/((R(L)^2)*(dt^2));
                        bc=0;
                        fc=1/dr^2+1/(R(L)*2*dr);
                        dc=1/dr^2-1/(R(L)*2*dr);
                        kc=0;
                        cc=0;
                        gc=0;
                        ac=0;

                    else  %internal nodes
                        ec=-2*( 1/dr^2 + 1/((R(L)^2)*(dt^2)) );
                        hc=1/((R(L)^2)*(dt^2));
                        bc=1/((R(L)^2)*(dt^2));
                        fc=1/dr^2+1/(R(L)*2*dr);
                        dc=1/dr^2-1/(R(L)*2*dr);
                        kc=0;
                        cc=0;
                        gc=0;
                        ac=0;
                    end

                    %small letters, elements of L and U

                  am(L)=ac;
                  bm(L)=(bc-psi*cc*fm(inw)-am(L)*fm(isw))/(1.d0-psi*fm(iw)*fm(inw));
                  cm(L)=cc-bm(L)*fm(iw);
                  dm(L)=(dc-am(L)*(2.*psi*gm(isw)+hm(isw))-bm(L)*gm(iw)) ...
                       / (1.+2.*psi*gm(is));
                  em(L)=ec+am(L)*(psi*gm(isw)-km(isw))-bm(L)*hm(iw)+ ...
                       cm(L)*(2.*psi*fm(inw)-gm(inw)+psi*km(inw))+ ...
                       dm(L)*(2.*psi*gm(is)-fm(is));
                  if (abs(em(L))<1e-10);break;end
                  em(L) = 1/em(L);  %! attention keeps the INVERSE of EPSILON!
                  fm(L)=( fc-bm(L)*km(iw)-cm(L)* ...
                       (hm(inw)+2.*psi*fm(inw)+2.*psi*km(inw))) * em(L);
                  gm(L)=(gc-dm(L)*hm(is)) * em(L);
                  hm(L)=(hc-dm(L)*(psi*gm(is)+km(is))) * em(L);
                  km(L)=kc * em(L);
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
                    
                    if  j==2 | j==Nr+1
                        L=ii(i,j);
                        
                        qrr(L) = u(L)-rhs(L);


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
        
            

            %Back-Front Substitution

            backfron();

            for L=1:ismax
                u(L)=u(L)+du(L);
            end
            
            
        end


    end

%% Back-Front Sustitution 
%A'Du = b - Au
%Calculation of du
 
    function backfron()
        
        for L=1:ismax+1
            aux(L) = 0;
            du(L)= 0;
        end

        %Front Substitution

        for i=2:Nt+1
            for j=2:Nr+1

                L=ii(i,j);
                iw=ii(i-1,j);
                is=ii(i,j-1);
                inw=ii(i-1,j+1);
                isw=ii(i-1,j-1);
                comsw=am(L)*aux(isw);
                comw=bm(L)*aux(iw);
                comnw=cm(L)*aux(inw);
                coms=dm(L)*aux(is);
                
                    
                aux(L)=(qrr(L)-comsw-comw-comnw-coms)*em(L);
  
            end
        end
        
        %Back Substitution

        for i=Nt+1:-1:2
            for j=Nr+1:-1:2
                L=ii(i,j);
                ie=ii(i+1,j);
                in=ii(i,j+1);
                ine=ii(i+1,j+1);
                ise=ii(i+1,j-1);
                comn=fm(L)*du(in);
                comse=gm(L)*du(ise);
                come=hm(L)*du(ie);
                comne=km(L)*du(ine);
                
                du(L)=aux(L)-comn-comse-come-comne;
                
            end
        end
    end

end