clear;
// Rad
m1= 0.07;
r = 0.04;
J1 = 1/2*m1*r*r;

// Aufbau
m2 = 1.5;
l =  0.03;
b = 0.07;
h = 0.1;
J2 = 1/12*m2*(b*b+h*h);
g=9.81;

P = 80 ;
Kd = 5;
Kom = 3;
Ki = 1;

phisoll=0.0;

function f=integrieren(t, y)
    
    phi1 = y(1);
    om1  = y(2);
    phi2 = y(3);
    om2  = y(4);
    eInteg=y(5);

    e      =  phisoll + phi2;       
    M = P * phi2 + Kd*om2 + Kom * om1 +Ki*eInteg;
    
    if(M>0.8)
        M=0.8;
    end
    if(M<-0.8)
        M=-0.8;
    end    
    A = l*sin(phi2);
    B = l*cos(phi2);

    MatrixA = [  (-m1*r*r)/J1 - 1, 0,                  1,              0;
                 0,               -1,                  0,              1;
                (-m2*r*r)/J1,     0, (-m2*B*B)/J2 - 1, (-m2*A*B)/J2;
                0,                0, (-m2*A*B)/J2,     (-m2*A*A)/J2 - 1  ];

    Matrixb = [m1*r*M/J1; 
                  -m1*g; 
               (m2*r*M)/J1 - (m2*B*M)/J2 - (m2*A*om2*om2);
              (-m2*A*M)/J2 + m2*B*om2*om2 - m2*g ];

    FZ = inv(MatrixA)*Matrixb;

    Fx1 = FZ(1);
    Fx2 = FZ(3);
    Fy2 = FZ(4);

    f(1) = om1;
    f(2) = (r*Fx1 + M)/J1;
    f(3) = om2;
    f(4) = (A*Fy2 + B*Fx2 - M)/J2;
    f(5) = e;
    
endfunction

t0 = 0;
y0 = [0,0,0.2,0,0]';
t = linspace(0,1,100);
y  = ode(y0,t0,t,integrieren);
plot(t,y(1,:)',t,y(3,:)');
legend("phi1","phi2",pos = 2) ;  

//Optimierung
function err=berechneFehler(x)
    P=x(1);
    Kd=x(2);
    Kom=x(3);
    Ki=x(4);
    t = linspace(0,1,100);
    y0 = [0,0,0.2,0,0]';
    t0 = 0;
    y  = ode(y0,t0,t,integrieren);

    err = sum(y(2,:).*y(2,:));
endfunction

P = 80 ;
Kd = 5;
Kom = 3;
Ki = 1;

x0=[P;Kd;Kom;Ki]; // Startpunkt eventuell zufaellig waehlen
ug = [0;0;0;0]; // untere Schranke fuer x
og = [1000;100;100;100]; // obere Schranke fuer x
[fopt, xopt] = optim(list(NDcost,berechneFehler),'b',ug,og,x0);

P=xopt(1);
Kd=xopt(2);
Kom=xopt(3);
Ki=xopt(4);

t0 = 0;
y0 = [0,0,0.2,0,0]';
t = linspace(0,1,100);
y  = ode(y0,t0,t,integrieren);
xset("window",2);
clf;
plot(t,y(1,:)',t,y(3,:)');
legend("phi1","phi2",pos = 2) ;

