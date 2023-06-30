accDir = pwd ;
addpath( genpath( [ accDir '/../source'] ) );
%Iteratively solves oscillator equation and non linear Vander Pol equation
loadParametersEx3
qsol = [q0;dq0];
Usol = [Z0;dZ0;X0;dX0];
% Time loop
n = 1; qtol = 1e-6; 
options = odeset('RelTol',1e-10);
itermax = 10; NumbIter = [];
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
            [t,Unp1] = ode45(@(t,Un) springBeam(t, Un,m,ku,cu,kq,cq,l,d,...
                rhoFluid,cL0,cD0,C2,va, qnp1k), t, [Zn; dZn; Xn; dXn; qnp1k], options) ;
            Udot = springBeam(t, Unp1(end, :),m,ku,cu,kq,cq,l,d...
                ,rhoFluid,cL0,cD0,C2,va, qnp1k);
            ddZnp1 = Udot(2);
            Usol(1:4, n+1) = Unp1(end, 1:4);
            % Solving qnp1 = ode(Ynp1, qn)
            %Using funcvanderpol
            [t,qnp1] = ode45(@(t, q) funcvanderpol(t,q, cq, kq, ddZnp1, d, A),t,[qn,dqn], options); %ode45(@(t,qn)funcvanderpol(t,qn,ddYn),t,[qn,dqn]); 
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

figure(1)  
hold on, grid on
yyaxis left
plot(times(1:spanPlot:end), Usol(1,1:spanPlot:end-1) ,'b-s', 'linewidth', lw,'markersize', ms+4)
xlabel('times')
ylabel('Ux')
hold on
yyaxis right
plot(times(1:spanPlot:end), qsol(1,1:spanPlot:end-1) ,'r-s', 'linewidth', lw,'markersize', ms+4 )
ylabel('q')

figure(2)
plot(times(1:spanPlot:end), Usol(3,1:spanPlot:end-1) ,'b-s', 'linewidth', lw,'markersize', ms+4)
xlabel('times')
ylabel('Ux')

% title(sprintf('m = %d, cu = %d, ku = %d, cq = %d, kq = %d', m, cu, ku, cq, kq))
%title(sprintf('Iterative Partitionned solution with m = %d, cu = %d, ku = %d, cq = %d, kq = %d', m, cu, ku, cq, kq))
%title(sprintf('Monolithic VS Iterative Partitionned solution'))
% legend('Y Monolithic  Newmark','Y Iterative Newmark','q Monolithic  Newmark','q Iterative Newmark');

