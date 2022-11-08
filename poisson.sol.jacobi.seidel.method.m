clear;clc
%% Mesh construction
IM=input('give number of nodes in x direction: ');
JM=input('give number of nodes in y direction: ');
Ntot=IM*JM;

len_x=input('give xmax: ');
len_y=input('give ymax: ');

dx=len_x/(IM-1); dy=len_y/(JM-1);

%% BOUNDARY CONDITIONS

    % Initiallization of u,A
    u=zeros(Ntot,1);
    A=zeros(Ntot,5);

    %right boundary
    m=0;
    for k=Ntot-JM+1:Ntot
        I=IM;
        u(k)=(dx*I)^2+(dy*m)^2;
        m=m+1;
        if m>JM-1 break;end
    end

    %left boundary
    m=0;
    for k=1:JM
       I=0;
       u(k)=(dx*I)^2+(dy*m)^2;
       m=m+1;
       if m>JM-1 break;end
    end

    %upper boundary
    m=0;
    for k=IM:JM:Ntot
       J=JM-1;
       u(k)=(dx*m)^2+(dy*J)^2;
       m=m+1;
       if m>IM-1 break;end
    end    

    %lower boundary
    m=0;
    for k=1:JM:Ntot-JM
       J=0;
       u(k)=(dx*m)^2+(dy*J)^2;
       m=m+1;
       if m>IM-1 break;end
    end  


%% Sparse matrix A (5 diagonals)
for i=2:IM-1
    for j=2:JM-1
    
    k=(i-1)*JM+j;
    A(k,1)=1/dy^2;
    A(k,2)=1/dx^2;
    A(k,3)=-2/dx^2-2/dy^2;
    A(k,4)=A(k,2);
    A(k,5)=A(k,1);
        
    end
end

%residual preallocation
r = zeros(Ntot,1);

q=4;
maxiter=1e6;

tol_1=1e-6; %tolerance for residual
tol_2=1e-10; %tolerance for unew-uold

omega=input('give relaxation omega: ');
iteration_method=input('Type 1 for Jacobi or 2 for Gauss Siedel method: ');

tic;

if iteration_method==1
    
    %% JACOBI
    res_ev=[];
    iterations=[];

    for p=1:maxiter+1
        uold=u;

        for i=2:IM-1
            for j=2:JM-1
                k=(i-1)*JM+j;

                %Delta formulation with relaxation
                du(k)= (1/A(k,3)) * ( q - A(k,1)*uold(k-JM) - A(k,2)*uold(k-1) - A(k,4)*uold(k+1) - A(k,5)*uold(k+JM) )-uold(k);
                u(k)=du(k)+uold(k);
                u(k)=omega*u(k)+(1-omega)*uold(k);

            end
        end

        for i=2:IM-1
            for j=2:JM-1
                k=(i-1)*JM+j;

                %Residual
                r(k)=A(k,1)*u(k-JM)+A(k,2)*u(k-1)+A(k,3)*u(k)+A(k,4)*u(k+1)+A(k,5)*u(k+JM)-q;

            end
        end

        res_ev(p)=norm(r);
        iterations(p)=p;

        err  = abs(norm(u-uold));
        res = err / (norm(u));

        if (norm(r)<tol_1)
            fprintf('Iteration stopped with tol_1 \n')
            fprintf('Total number of iterations: %d \n',p)
            break;end
        if (res < tol_2)
            fprintf('Iteration stopped with tol_2 \n')
            fprintf('Total number of iterations: %d \n',p)
            break;end

    end

elseif iteration_method==2
    
    %% GAUSS SEIDEL
    res_ev=[];
    iterations=[];

    for p=1:maxiter+1

        uold=u;

        for i=2:IM-1
            for j=2:JM-1
                k=(i-1)*JM+j;

                %Delta formulation with relaxation
                du(k)= (1/A(k,3)) * ( q - A(k,1)*u(k-JM) - A(k,2)*u(k-1) - A(k,4)*uold(k+1) - A(k,5)*uold(k+JM) )-uold(k);
                u(k)=du(k)+uold(k);
                u(k)=omega*u(k)+(1-omega)*uold(k);
            end
        end

        for i=2:IM-1
            for j=2:JM-1
                k=(i-1)*JM+j;

                %Residual
                r(k)=A(k,1)*u(k-JM)+A(k,2)*u(k-1)+A(k,3)*u(k)+A(k,4)*u(k+1)+A(k,5)*u(k+JM)-q;

            end
        end


        res_ev(p)=norm(r);
        iterations(p)=p;


        res  = norm(u-uold)/norm(u);

       
        if (norm(r)<tol_1)
            fprintf('Iteration stopped with tol_1 \n')
            fprintf('Total number of iterations: %d \n',p)
            break;end
        if (res < tol_2)
            fprintf('Iteration stopped with tol_2 \n')
            fprintf('Total number of iterations: %d \n',p)
            break;end


    end
    
end
toc;

%% Residual drop factor
res_drop=-log10(res_ev(end)/res_ev(1));
fprintf('Residual drop factor(first to last iteration): %.3g \n %f \n',res_drop)

%% Grid coordinates
xcoor=0:dx:len_x;
ycoor=0:dy:len_y;
[X,Y]=meshgrid(xcoor,ycoor);

figure(1)
Z=zeros(IM);
mesh(X,Y,Z)
xlabel('x')
ylabel('y')
zlabel('z')
title("Mesh")

%% Analytical sol
for i=1:IM
    for j=1:JM
        u_exact(i,j)=X(i,j)^2+Y(i,j)^2;
    end
end

figure(1)
contour(X,Y,u_exact,'ShowText','on')

%% Results figures

n=JM;
w=1;
k=1;
 for j=1:JM
    uplot(1:JM,k)=u(w:n,1);
   
    n=n+JM;
    k=k+1;
    w=w+JM;
    
 end

figure(2)
subplot(1,2,1)
surf(X,Y,uplot)
xlabel('x')
ylabel('y')
zlabel('z')
title("3D Solution")

subplot(1,2,2)
contour(X,Y,uplot,'ShowText','on')
xlabel('x')
ylabel('y')
title("2D Solution")

figure(3)
subplot(1,2,1)
plot(res_ev,iterations)
ylabel('residual')
xlabel('#iterations')
title("Residual-Number of iterations")
hold on

subplot(1,2,2)
semilogy(res_ev,iterations)
ylabel('residual')
xlabel('#iterations')
title("Residual-Number of iterations (logarithmic scale)")
hold off

figure(4)
contour(X,Y,u_exact,"r",'ShowText','on')
hold on
contour(X,Y,uplot,"g",'ShowText','on')
xlabel('x')
ylabel('y')
legend('Analytical','Iterative method')
hold off
