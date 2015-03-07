A=[0 0];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5);
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[0 sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5);
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[0 -sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[0 2*sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[0 -2*sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

% right
A=[1.5 sqrt(3)/2];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[1.5 -sqrt(3)/2];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[1.5 1.5*sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[1.5 -1.5*sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);



A=[3 0];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[3 sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[3 -sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);



% left
A=[-1.5 sqrt(3)/2];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[-1.5 -sqrt(3)/2];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[-1.5 1.5*sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[-1.5 -1.5*sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);



A=[-3 0];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[-3 sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

A=[-3 -sqrt(3)];
N=1;
B=[];C=[];
B=[A(1)+N A(1)+N/2 A(1)-N/2 A(1)-N A(1)-N/2 A(1)+N/2 A(1)+N];
C=[A(2) A(2)+N*sqrt(3)/2 A(2)+N*sqrt(3)/2 A(2) A(2)-N*sqrt(3)/2 A(2)-N*sqrt(3)/2 A(2)];
plot(B,C,'k','LineWidth', 1.5)
hold on
plot(A(1),A(2),'ro','LineWidth', 1.5)
users = brownian(1, 5, A(1)+A(2)*1j, 1);
plot(users, 'bx', 'LineWidth', 1);

axis equal

