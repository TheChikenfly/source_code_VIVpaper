function IterativeODE3test()
%Iteratively solves beam subject to uniform lift force and non linear Vander Pol equation
% flow parameters
%close all;
booldrag = true;
global qvect;
% time resolution
dt = 0.00028; 
finalTime =0.2;
A = 12; epsilon = 0.3; St = 0.16;
D = 0.01; l = 0.2;
%D = 1e-1;
CL0 = 0.3; CD0 = 0;
if booldrag
    CD0 = 1.2; % CAN BE CHANGED TO 0 
end
va = 0.4; rhoFluid = 1020;
%Beam parameters
rho = rhoFluid ; %Neutraly buoyant
E = 5e8; I = pi*D^4/64;
% 1D Oscillator parameters
m = rho*pi*(D^2)*l/(4*2); ku = 3*E*I/(l^3); cu = 0; %fu = 0;
B = 2*pi*St*va/D;
cq = epsilon*B; % Vpr = Va here since rod has 2D motion on plane perp to flow
kq = B^2; 
C2 = A/D; 
%Initial conditions
q0 = 2; dq0 = 0; ddq0 = 0; 
Z0 = 0; dZ0 = 0; ddZ0 = 0;%0.5 * rhoFluid * CL0 * D * l * 3 * va^2 * q0/(m*8*2); %0; First acceleration
X0 = 0; dX0 = 0; ddX0 = 0;%0.5 * rhoFluid * CD0 * d * l * 3 * va^2 * q0/(8*2);
qsol = [q0;dq0];
%Usol = [Z0;dZ0;ddZ0 ;X0;dX0;ddX0];
Usol = [Z0;dZ0;X0;dX0];
% Time loop
n = 1; qtol = 1e-6; 
options = odeset('RelTol',1e-10);
itermax = 50; NumbIter = [];
while n*dt <= finalTime 
    N= 5; h= dt/N; tn = n*dt;% Time step
    t = tn:h:tn + dt;
    qn = qsol(1, n); dqn = qsol(2, n);
    Zn = Usol(1, n); dZn = Usol(2, n); %ddZn = Usol(3, n);
    Xn = Usol(3, n); dXn = Usol(4, n); %ddXn = Usol(6, n);
    qnp1k = qsol(1, n); % candidate for time n+1
    flag = true;
    k = 0;
        while flag
            % Solving displacement Ynp1 = ode(Yn, qn)
            [t,Unp1] = ode45(@(t,Un) springBeam(t, Un,m,ku,cu,kq,cq,l,D,...
                rhoFluid,CL0,CD0,C2,va, qnp1k), t, [Zn; dZn; Xn; dXn; qnp1k], options) ;
            Udot = springBeam(t, Unp1(end, :),m,ku,cu,kq,cq,l,D...
                ,rhoFluid,CL0,CD0,C2,va, qnp1k);
            ddZnp1 = Udot(2);
            Usol(1:4, n+1) = Unp1(end, 1:4);
            % Solving qnp1 = ode(Ynp1, qn)
            %Using funcvanderpol
            [t,qnp1] = ode45(@(t, q) funcvanderpol(t,q, cq, kq, ddZnp1, D, A),t,[qn,dqn], options); %ode45(@(t,qn)funcvanderpol(t,qn,ddYn),t,[qn,dqn]); 
            qsol(1:2, n+1) = qnp1(end, 1:2); % The N+1 component
            qnp1kp1 = qsol(1, n+1);
            if abs((qnp1kp1 - qnp1k )/qnp1k) < qtol 
                flag = false;
                qsol(1:2, n+1) = qnp1(end, 1:2); % The N+1 component
            elseif k == itermax
                fprintf('maxiter reached')
                stop;
            end
            qnp1k = qnp1kp1; 
%             figure(8)
%             plot(n, qnp1kp1 ,'k');
%             hold on
            k = k+1; 
        end
    NumbIter(n) = k-1 ;
    n = n+1 
end

% T = table(effectiveload(1:2,:),dUt(1:2, :))
% writetable(T,'NewmarkValues.dat', 'WriteRowNames',true)  
% type 'NewmarkValues.dat';
N = finalTime/dt; % nombre de valeurs à afficher %all values
times = 0:dt:N*dt-dt;%finalTime-dt;
spanPlot = 20 ; lw = 2.0 ; ms = 11 ; plotfontsize = 20 ;
figure(10)  
hold on, grid on
yyaxis left
length(times(1:spanPlot:end))
length(Usol(1,1:spanPlot:end-1))
plot(times(1:spanPlot:end), Usol(1,1:spanPlot:end-1) ,'b-o', 'linewidth', lw,'markersize', ms )
xlabel('times')
ylabel('Uz')
hold on
yyaxis right
plot(times(1:spanPlot:end), qsol(1,1:spanPlot:end-1) ,'r-o', 'linewidth', lw,'markersize', ms )
% axis auto
ylabel('q')
title(sprintf('Iterative Partitionned solution with m = %d, cu = %d, ku = %d, cq = %d, kq = %d', m, cu, ku, cq, kq))
%legend('Y Monolithic  Newmark','Y Iterative Newmark','q Monolithic  Newmark','q Iterative Newmark');

figure(2) 
plot(times(1:spanPlot:end), Usol(3,1:spanPlot:end-1) ,'b-o', 'linewidth', lw,'markersize', ms )

figure(4)
plot(NumbIter)
xlabel('n')
ylabel('k')
title('Number of iteration before convergence')
end


