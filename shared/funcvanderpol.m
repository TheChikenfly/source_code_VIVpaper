function dq=funcvanderpol(t,qn, c, k, ddYnp1, D, A)
    qdot = qn(2);
    qdotdot = c*(1-qn(1)^2)*qn(2)-k*qn(1)+A*ddYnp1/D;
    dq = [qdot; qdotdot]; 
end

