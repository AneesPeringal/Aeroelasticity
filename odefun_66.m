function v=odefun_66(t,x,A,Fq1,Fq2)
v=A*x+[0 ;0; 0; Fq1  ].*x.^3+[0 ;0 ;0 ;Fq2].* (x(1)-x(2)-x(3))^3;
end