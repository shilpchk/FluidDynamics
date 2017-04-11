function sbnewt
%Changed l2 to match snewt

global C del L1 L2 R1 R2 B d U eps amp l2 sul ss

%   Newtonian integration

deli=0.2;
dele=0.1;
nd=20;
ddel=(deli-dele)/nd;
del=deli;
C=1;

lamL=(roots([1 0 -1/C (2-del-del^2)/C]));
[dl,I]=sort(real(lamL));
L1=lamL(I(2));
L2=lamL(I(3));
%[L1 L2]
lamR=(roots([1 0 -1/C (2*del^2-del-1)/del^3/C]));
[dl,I]=sort(real(lamR));
R1=lamR(I(1));
R2=lamR(I(2));
%[R1 R2]

%options=bvpset
options = bvpset('Vectorized','on','Nmax',4000,'FJacobian',...
     @nogJacN,'BCJacobian',@nogBCJacN); %,'AbsTol',1e-8,'RelTol',1e-5);

%Calculate initial Newtonian film shape
%Strangely enough this works better than in newt2, allows greater $x$
%range
solinit = bvpinit(linspace(-7.5,4.5,200),@noginitN);
sol = bvp4c(@nogodeN,@nogbcN,solinit,options);

figure(1)
hold on, plot(sol.x,sol.y(1,:)), hold off

%Now reduce del for Newtonian film
for n=1:nd
del = del-ddel;
lamL=(roots([1 0 -1/C (2-del-del^2)/C]));
[dl,I]=sort(real(lamL));
L1=lamL(I(2));
L2=lamL(I(3));
lamR=(roots([1 0 -1/C (2*del^2-del-1)/del^3/C]));
[dl,I]=sort(real(lamR));
R1=lamR(I(1));
R2=lamR(I(2));
[L1 L2 R1 R2]
 
sol = bvp4c(@nogodeN,@nogbcN,sol,options);
%Plot n film shapes
  if (mod(n,10)==0)
    figure(1)
    hold on,plot(sol.x,sol.y(1,:))
    hold off
  end
title('Newtonian solutions for decreasing delta')
end

%   non-Newtonian integration

B=1e-3;
eps=1e-6;

Uflux
[B U d]
h=sol.y(1,:);

%Don't forget the Newtonian solution
hN=h;
xN=sol.x;

[Gs,F]=Gamit(h);
Y=h-B./Gs;

%Plot final h again for Newtonian case. very sensible.
%figure(7)
%hold on, plot(sol.x,h), hold off
%axis([min(sol.x) max(sol.x) 0 max(h)])

p=polyB(1);
lamL=(roots(p));
[dl,I]=sort(real(lamL));
L1=lamL(I(2));
L2=lamL(I(3));
%[L1 L2]

p=polyB(del);
lamR=(roots(p));
[dl,I]=sort(real(lamR));
R1=lamR(I(1));
R2=lamR(I(2));
%[R1 R2]

%options=bvpset;
options = bvpset('Vectorized','on','Nmax',4000,'FJacobian',...
     @nogJac,'BCJacobian',@nogBCJac,'AbsTol',1e-6,'RelTol',1e-4);

sol = bvp4c(@nogode,@nogbc,sol,options);

h=sol.y(1,:);
[Gs,F]=Gamit(h);
Y=h-B./Gs;

%Plot the first non-Newtonian solution and compare with Newtonian
figure(2)
hold on,plot(sol.x,h,sol.x,Y, xN,hN)
hold off
title('Comparison of Newtonian, first non-Newtonian, B=1e-3, and yield surface')

options = bvpset('Vectorized','on','Nmax',8000,'FJacobian',@snogJac); %,...
%'AbsTol',1e-6,'RelTol',1e-4);
%,'BCJacobian',@snogBCJac);

l2=0.001^2;
 
U=(1-del^3)/3/(1-del);
d=1/3-U;
amp=1e-6;

sul=sol;

sol.y(1,:)=sol.y(2,:);
sol.y(2,:)=sol.y(3,:);
%Approximating the derivatives
sol.y(3,:)=diff([sol.y(3,:) 0])./diff([sol.x sol.x(end)*2-sol.x(end-1)]);
hppp=diff([sol.y(3,:) 0])./diff([sol.x sol.x(end)*2-sol.x(end-1)]);

H=sul.y(1,:);
hp=sol.y(2,:);
[Gs,F,Fh,FG,Gh]=Gamit(H);
Q=Fh.*(d+U*H)./F+(C*hppp-(1+C*l2)*hp).*(F+FG.*abs((d+U*H)./F));
sol.y(4,:)=Q;
sol.y(1:4,:)=sol.y(1:4,:)/sol.y(3,1)*amp;
%sol.y(3.1)
sol.parameters=0.001;
 
' Stability integration'
 
sol = bvp4c(@snogode,@snogbc,sol,options);
 
ss(1)=sol.parameters; sol.parameters
 
figure(3)
plot(sol.x,sol.y(1,:)/max(abs(sol.y(1,:))))

%Here I've changed l2 (was 0.015*n), made changes to figure(3) below
for n=1:50
l2=(0.01*n)^2;
amp=1e-5/(10+n+(n/30)^4);
sol = bvp4c(@snogode,@snogbc,sol,options);
ss(n+1)=sol.parameters; sol.parameters
figure(3)
hold on,plot(sol.x,sol.y(1,:)/max(abs(sol.y(1,:)))),hold off
%hold on,plot(sol.x,sol.y(4,:)/max(abs(sol.y(4,:)))),hold off
end
 
%Was 0.015*[ 
figure(4)
plot(0.0125*[0:50],ss)

B=0;
l=0;
for n=1:15
B=B+0.75e-3;

Uflux

p=polyB(1);
lamL=(roots(p));
[dl,I]=sort(real(lamL));
L1=lamL(I(2));
L2=lamL(I(3));
%[L1 L2]

p=polyB(del);
lamR=(roots(p));
[dl,I]=sort(real(lamR));
R1=lamR(I(1));
R2=lamR(I(2));
%[R1 R2]

sol=sul;

options = bvpset('Vectorized','on','Nmax',4000,'FJacobian',...
     @nogJac,'BCJacobian',@nogBCJac); %,'AbsTol',1e-6,'RelTol',1e-4);
sol = bvp4c(@nogode,@nogbc,sol,options);

h=sol.y(1,:);
[Gs,F]=Gamit(h);
Y=h-B./Gs;

[B U d]
% I'm glad this is commented out
%l=l+1;
%if l==2,
%  l=0;
%Plot the non-Newtonian solutions, height and yield surface
if (mod(n,5)==0)|(n==1)
  figure(5)
  hold on,plot(sol.x,h,sol.x,Y)
  hold off
end
title('Plot of height and yield surface for increasing B')
%end

options = bvpset('Vectorized','on','Nmax',8000,'FJacobian',@snogJac);
%,'BCJacobian',@snogBCJac);

l2=0.001^2;
 
U=(1-del^3)/3/(1-del);
d=1/3-U;
amp=1e-6;

sul=sol;

sol.y(1,:)=sol.y(2,:);
sol.y(2,:)=sol.y(3,:);
sol.y(3,:)=diff([sol.y(3,:) 0])./diff([sol.x sol.x(end)*2-sol.x(end-1)]);
hppp=diff([sol.y(3,:) 0])./diff([sol.x sol.x(end)*2-sol.x(end-1)]);
H=sul.y(1,:);
hp=sol.y(2,:);
[Gs,F,Fh,FG,Gh]=Gamit(H);
Q=Fh.*(d+U*H)./F+(C*hppp-(1+C*l2)*hp).*(F+FG.*abs((d+U*H)./F));
sol.y(4,:)=Q;
sol.y(1:4,:)=sol.y(1:4,:)/sol.y(3,1)*amp;
%sol.y(3.1)
sol.parameters=0.001;
 
' Stability integration'
 
sol = bvp4c(@snogode,@snogbc,sol,options);
 
ss(1)=sol.parameters; sol.parameters
 
figure(3)
plot(sol.x,sol.y(1,:)/max(abs(sol.y(1,:))))

for n=1:50
l2=(0.01*n)^2;
amp=1e-5/(10+n+(n/30)^4);
sol = bvp4c(@snogode,@snogbc,sol,options);
ss(n+1)=sol.parameters; sol.parameters
figure(3)
hold on,plot(sol.x,sol.y(1,:)/max(abs(sol.y(1,:)))),hold off
end
 

if (mod(n,5)==0)|(n==1)
  figure(4)
  hold on,plot(0.0125*[0:50],ss)
  hold off
end
title('The useful solution')

end
 
function res = snogbc(ya,yb,lam)
global C del B d U l2 d U  L1 L2 R1 R2  amp sul

%bcL=L1*L2*(ya(1)-1)+ya(3)-(L1+L2)*ya(2);
%bcR=R1*R2*(yb(1)-del)+yb(3)-(R1+R2)*yb(2);

H=1;
[Gs,F,Fh,FG,Gh]=Gamit(H);
c(1)=U*H+d+FG; c(2)=0; c(3)=-F*(1+2*C*l2)-FG*(1+C*l2);
c(4)=Fh-U; c(5)=lam+l2*F*(1+C*l2);
L=roots(c);
[dl,I]=sort(real(L));
cL1=L(I(3));
cL2=L(I(4));
ya7=((1+C*l2)*ya(2)+(ya(4)-Fh*(d+U*H)/F*ya(1))/(F+Gs*FG))/C;
 
H=del;
[Gs,F,Fh,FG,Gh]=Gamit(H);
c(1)=U*H+d+FG; c(2)=0; c(3)=-F*(1+2*C*l2)-FG*(1+C*l2);
c(4)=Fh-U; c(5)=lam+l2*F*(1+C*l2);
L=roots(c);
[dl,I]=sort(-real(L));
cR1=L(I(3));
cR2=L(I(4));
yb7=((1+C*l2)*yb(2)+(yb(4)-Fh*(d+U*H)/F*yb(1))/(F+Gs*FG))/C;
                                                                              
%[cL1 cL2 cR1 cR2]
                                                                             
bc1L=cL1*cL2*ya(1)+ya(3)-(cL1+cL2)*ya(2);
bc1R=cR1*cR2*yb(1)+yb(3)-(cR1+cR2)*yb(2);
 
bc2L=cL1*cL2*ya(2)+ya7-(cL1+cL2)*ya(3);
bc2R=cR1*cR2*yb(2)+yb7-(cR1+cR2)*yb(3);
 
res = [  ya(3)-amp
         ya(1)
         ya(2)
         yb(1)
         yb(2)];
res = [  ya(3)-amp
         bc1L
         bc2L
         bc1R
         bc2R];


function dydx = snogode(x,y,lam)
global C del B d U l2 d U sul

yi=deval(sul,x);

H=yi(1,:); Hp=yi(2,:); Hpp=yi(3,:);
h=y(1,:); hp=y(2,:); hpp=y(3,:); Q=y(4,:);

[Gs,F,Fh,FG,Gh]=Gamit(H);

dydx = [          y(2,:)
                  y(3,:)
  ((1+C*l2)*hp+(Q-Fh.*(d+U*H).*h./F)./(F+FG.*Gs))/C
     U*hp-lam*h+l2*F.*(C*hpp-(1+C*l2)*h)];

 
function [jac,dfdp] = snogJac(x,y,lam)
global C del l2 d U sul
 
yi=deval(sul,x);
 
H=yi(1,:); Hp=yi(2,:); Hpp=yi(3,:);
 
h=y(1,:); hp=y(2,:); hpp=y(3,:); hppp=y(4,:);

[Gs,F,Fh,FG,Gh]=Gamit(H);
 
jac = [   0 1 0 0
          0 0 1 0
 -Fh.*(d+U*H)./F./(F+FG.*Gs)/C (1+C*l2)/C 0 1./(F+FG.*Gs)/C
  -lam-l2*F*(1+C*l2) U l2*F 0];

dfdp=[0
      0
      0
     -h];

function Uflux
global C del L1 L2 R1 R2 B d U eps

H=1;
T=4*eps*B; Z=B-eps; Y=H-(B+eps);
th1=asinh(-Z/sqrt(T)); th2=asinh((H-Z)/sqrt(T));
F1=((H-Z)^2+T)^(3/2)/6-(Z^2+T)^(3/2)/6-H^3/12+Y*H^2/4+...
     T/8*Z*(sinh(th2*2)-sinh(th1*2)+2*(th2-th1));
H=del;
T=4*eps*B; Z=B-eps; Y=H-(B+eps);
th1=asinh(-Z/sqrt(T)); th2=asinh((H-Z)/sqrt(T));
F2=((H-Z)^2+T)^(3/2)/6-(Z^2+T)^(3/2)/6-H^3/12+Y*H^2/4+...
     T/8*Z*(sinh(th2*2)-sinh(th1*2)+2*(th2-th1));

U=(F1-F2)/(1-del);
d=F1-U;

function p=polyB(H)
global C del L1 L2 R1 R2 B d U eps

T=4*eps*B; Z=B-eps; Y=H-(B+eps);
th1=asinh(-Z/sqrt(T)); th2=asinh((H-Z)/sqrt(T));
F=((H-Z)^2+T)^(3/2)/6-(Z^2+T)^(3/2)/6-H^3/12+Y*H^2/4+...
     T/8*Z*(sinh(th2*2)-sinh(th1*2)+2*(th2-th1));
dZ=-Z; dT=-2*T; dY=B+eps;
dth1=-1/cosh(th1)/sqrt(T)*dZ-tanh(th1)/2/T*dT;
dth2=-1/cosh(th2)/sqrt(T)*dZ-tanh(th2)/2/T*dT;
dFG=T/4*Z*(cosh(th2*2)*dth2-cosh(th1*2)*dth1+dth2-dth1)+...
    dT/8*Z*(sinh(th2*2)-sinh(th1*2)+2*(th2-th1))+...
    T/8*dZ*(sinh(th2*2)-sinh(th1*2)+2*(th2-th1))+...
 ((H-Z)^2+T)^(1/2)/4*(dT-dZ*2*(H-Z))-(Z^2+T)^(1/2)/4*(2*Z*dZ+dT)+dY*H^2/4;

dth2=1/cosh(th2)/sqrt(T);
dFH=T/4*Z*(cosh(th2*2)+1)*dth2+((H-Z)^2+T)^(1/2)/2*(H-Z)+Y*H/2;

p(1)=(1+(d+H*U)/F^2*dFG)*C;
p(3)=-p(1)/C;
p(4)=(d+U*H)*dFH/F^2-U/F;

function [Gs,F,Fh,FG,Gh]=Gamit(h)
global C del L1 L2 R1 R2 B d U eps

Gs=3./h.^3.*(d+U*h);
for n=1:20
T=4*eps*B./Gs.^2;
Z=(B-eps)./Gs;
Y=h-(B+eps)./Gs;
th1=asinh(-Z./sqrt(T));
th2=asinh((h-Z)./sqrt(T));
F=((h-Z).^2+T).^(3/2)/6-(Z.^2+T).^(3/2)/6-h.^3/12+Y.*h.^2/4+...
     T/8.*Z.*(sinh(th2*2)-sinh(th1*2)+2*(th2-th1));
dZ=-Z./Gs;
dT=-2*T./Gs;
dY=(B+eps)./Gs.^2;
dth1=-1./cosh(th1)./sqrt(T).*dZ-tanh(th1)/2./T.*dT;
dth2=-1./cosh(th2)./sqrt(T).*dZ-tanh(th2)/2./T.*dT;
dF=T/4.*Z.*(cosh(th2*2).*dth2-cosh(th1*2).*dth1+dth2-dth1)+...
    dT/8.*Z.*(sinh(th2*2)-sinh(th1*2)+2*(th2-th1))+...
    T/8.*dZ.*(sinh(th2*2)-sinh(th1*2)+2*(th2-th1))+...
     ((h-Z).^2+T).^(1/2)/4.*(dT-dZ*2.*(h-Z))-...
         (Z.^2+T).^(1/2)/4.*(2*Z.*dZ+dT)+dY.*h.^2/4;
funk=abs(d+h*U)-abs(F).*Gs;
dfunk=-abs(F)-Gs.*dF.*sign(F);
dG=-funk./dfunk; 
%if n==21,
%  dG=dG/2;
%end
Gs=Gs+dG;
%[mean(abs(dG))]
if mean(abs(dG))<1e-12,
break
end
end

if mean(abs(dG))>1e-6,
mean(abs(dG))
end

dth2=1./sqrt(T)./cosh(th2);
Fh=(h-Z).*((h-Z).^2+T).^(1/2)/2+Y.*h/2+T/4.*Z.*(cosh(th2*2)+1).*dth2;
funkh=U*sign(d+h*U)-Gs.*Fh.*sign(F);
Gh=-funkh./dfunk;
FG=dF;

function dydx = nogode(x,y)
global C del B d U 

h=y(1,:);
[Gs,F,Fh,FG,Gh]=Gamit(h);
Y=h-B./Gs;

dydx = [          y(2,:)
                  y(3,:)
           ((d+h*U)./F-1+y(2,:))/C];

function jac = nogJac(x,y)
global del C B d U

h=y(1,:);
[Gs,F,Fh,FG,Gh]=Gamit(h);

jac = [   0  1  0
          0  0  1
 (U./F-(d+h*U)./F.^2.*(Fh+FG.*Gh))/C 1/C 0];

% --------------------------------------------------------------------------
 
function res = nogbc(ya,yb)
global del L1 L2 R1 R2

bcL=L1*L2*(ya(1)-1)+ya(3)-(L1+L2)*ya(2);
bcR=R1*R2*(yb(1)-del)+yb(3)-(R1+R2)*yb(2);

res = [  ya(1)-1
         bcL
         bcR ];
 
function yinit = noginit(x)
global del
     yinit = [   1-(1+tanh(4*(x-2)))/2*(1-del)
     1./cosh(x).^2*(1-del)/2*4
                -sinh(x)./cosh(x).^3*(1-del)*16];
 
% --------------------------------------------------------------------------
 
function [dBCdya,dBCdyb] = nogBCJac(ya,yb)

global del L1 L2 R1 R2

dBCdya = [ 1  0  0
      L1*L2 -(L1+L2) 1
           0  0  0]; 
dBCdyb = [ 0  0  0
           0  0  0
      R1*R2 -(R1+R2) 1]; 

function dydx = nogodeN(x,y)
global C del

dydx = [              y(2,:)
                      y(3,:)
   (((1+del+del^2)*y(1,:)-del-del^2)./y(1,:).^3-1+y(2,:))/C];


% --------------------------------------------------------------------------
 
function res = nogbcN(ya,yb)
global del L1 L2 R1 R2

bcL=L1*L2*(ya(1)-1)+ya(3)-(L1+L2)*ya(2);
bcR=R1*R2*(yb(1)-del)+yb(3)-(R1+R2)*yb(2);

res = [  ya(1)-1
         bcL
         bcR ];
 
function yinit = noginitN(x)
global del
     yinit = [   1-(1+tanh(4*(x-2)))/2*(1-del)
     1./cosh(x).^2*(1-del)/2*4
                -sinh(x)./cosh(x).^3*(1-del)*16];

function jac = nogJacN(x,y)
global del C
jac = [   0  1  0
          0  0  1
 (-2*(1+del+del^2)./y(1,:).^3+3*del*(1+del)./y(1,:).^4)/C 1/C 0];
 
% --------------------------------------------------------------------------
 
function [dBCdya,dBCdyb] = nogBCJacN(ya,yb)

global del L1 L2 R1 R2

dBCdya = [ 1  0  0
      L1*L2 -(L1+L2) 1
           0  0  0]; 
dBCdyb = [ 0  0  0
           0  0  0
      R1*R2 -(R1+R2) 1]; 

