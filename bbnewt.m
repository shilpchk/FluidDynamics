function bbnewt
global C del L1 L2 R1 R2 B d U eps

%   Newtonian integration

del=0.2;
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

%options=bvpsetc
options = bvpsetc('Vectorized','on','Nmax',4000,'FJacobian',...
     @nogJacN,'BCJacobian',@nogBCJacN,'AbsTol',1e-8,'RelTol',1e-5);

solinit = bvpinit(linspace(-7,6,100),@noginitN);
sol = bvp4cc(@nogodeN,@nogbcN,solinit,options);

%   non-Newtonian integration

B=0.1;
eps=0.001;

Uflux
[B U d]
h=sol.y(1,:);
[Gs,F]=Gamit(h);
Y=h-B./Gs;

figure(1)
plot(sol.x,h)
axis([min(sol.x) max(sol.x) 0 max(h)])

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

%options=bvpsetc;
options = bvpsetc('Vectorized','on','Nmax',4000,'FJacobian',...
     @nogJac,'BCJacobian',@nogBCJac,'AbsTol',1e-6,'RelTol',1e-4);

sol = bvp4cc(@nogode,@nogbc,sol,options);

h=sol.y(1,:);
[Gs,F]=Gamit(h);
Y=h-B./Gs;

figure(1)
hold on,plot(sol.x,h,sol.x,Y)
hold off

l=0;
for n=1:9
B=B+0.1;

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

sol = bvp4cc(@nogode,@nogbc,sol,options);

h=sol.y(1,:);
[Gs,F]=Gamit(h);
Y=h-B./Gs;

[B U d]
%l=l+1;
%if l==2,
%  l=0;
figure(1)
hold on,plot(sol.x,h,sol.x,Y)
hold off
end
%end


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
[Gs,F]=Gamit(h);
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

res = [  ya(1)-1+1e-4
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

function sol = bvp4cc(ode, bc, solinit, options, varargin);
%BVP4C  Solve two-point boundary value problems for ODEs by collocation.
%   SOL = BVP4C(ODEFUN,BCFUN,SOLINIT) integrates a system of ordinary
%   differential equations of the form y' = f(x,y) on the interval [a,b],
%   subject to general two-point boundary conditions of the form
%   bc(y(a),y(b)) = 0. ODEFUN is a function of two arguments: a scalar X
%   and a vector Y. ODEFUN(X,Y) must return a column vector representing
%   f(x,y). BCFUN is a function of two vector arguments. BCFUN(YA,YB) must
%   return a column vector representing bc(y(a),y(b)). SOLINIT is a structure
%   with fields named
%       x -- ordered nodes of the initial mesh with
%            SOLINIT.x(1) = a, SOLINIT.x(end) = b
%       y -- initial guess for the solution with SOLINIT.y(:,i)
%            a guess for y(x(i)), the solution at the node SOLINIT.x(i)
%
%   BVP4C produces a solution that is continuous on [a,b] and has a
%   continuous first derivative there. The solution is evaluated at points
%   XINT using the output SOL of BVP4C and the function DEVAL:
%   YINT = DEVAL(SOL,XINT). The output SOL is a structure with
%       SOL.x  -- mesh selected by BVP4C
%       SOL.y  -- approximation to y(x) at the mesh points of SOL.x
%       SOL.yp -- approximation to y'(x) at the mesh points of SOL.x
%       SOL.solver -- 'bvp4c'
%
%   SOL = BVP4C(ODEFUN,BCFUN,SOLINIT,OPTIONS) solves as above with default
%   parameters replaced by values in OPTIONS, a structure created with the
%   BVPSET function. To reduce the run time greatly, use OPTIONS to supply
%   a function for evaluating the Jacobian and/or vectorize ODEFUN.
%   See BVPSET for details and SHOCKBVP for an example that does both.
%
%   SOL = BVP4C(ODEFUN,BCFUN,SOLINIT,OPTIONS,P1,P2...) passes constant, known
%   parameters P1, P2... to the functions ODEFUN and BCFUN, and to all
%   functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
%   no options are set.
%
%   Some boundary value problems involve a vector of unknown parameters p
%   that must be computed along with y(x):
%       y' = f(x,y,p)
%       0  = bc(y(a),y(b),p)
%   For such problems the field SOLINIT.parameters is used to provide a guess
%   for the unknown parameters. On output the parameters found are returned
%   in the field SOL.parameters. The solution SOL of a problem with one set
%   of parameter values can be used as SOLINIT for another set. Difficult BVPs
%   may be solved by continuation: start with parameter values for which you can
%   get a solution, and use it as a guess for the solution of a problem with
%   parameters closer to the ones you want. Repeat until you solve the BVP
%   for the parameters you want.
%
%   The function BVPINIT forms the guess structure in the most common
%   situations:  SOLINIT = BVPINIT(X,YINIT) forms the guess for an initial mesh X
%   as described for SOLINIT.x and YINIT either a constant vector guess for the %   solution or a function that evaluates the guess for the solution
%   at any point in [a,b]. If there are unknown parameters,
%   SOLINIT = BVPINIT(X,YINIT,PARAMS) forms the guess with the vector PARAMS of %   guesses for the unknown parameters.
%
%   Example
%         solinit = bvpinit([0 1 2 3 4],[1 0]);
%         sol = bvp4c(@twoode,@twobc,solinit);
%     solve a BVP on the interval [0,4] with differential equations and
%     boundary conditions computed by functions twoode and twobc, respectively. %     This example uses [0 1 2 3 4] as an initial mesh, and [1 0] as an initial %     approximation of the solution components at the mesh points.
%         xint = linspace(0,4);
%         yint = deval(sol,xint);
%     evaluate the solution at 100 equally spaced points in [0 4]. The first
%     component of the solution is then plotted with
%         plot(xint,yint(1,:));
%   For more examples see TWOBVP, MAT4BVP and SHOCKBVP.
%
%   See also BVPSET, BVPGET, BVPINIT, DEVAL, @.

%   BVP4C is a finite difference code that implements the 3-stage Lobatto
%   IIIa formula. This is a collocation formula and the collocation
%   polynomial provides a C1-continuous solution that is fourth order
%   accurate uniformly in [a,b]. Mesh selection and error control are based
%   on the residual of the continuous solution. Analytical condensation is
%   used when the system of algebraic equations is formed.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2001 The MathWorks, Inc.
%   $Revision: 1.13 $  $Date: 2001/03/01 14:46:57 $

true = 1;
false = ~true;

% check input parameters
if nargin <3
  error('Not enough input arguments')
end

if ~isstruct(solinit)
  error('The initial profile must be provided as a structure')
elseif ~isfield(solinit,'x')
  error(['The field ''x'' not present in ''' inputname(3) ''''])
elseif ~isfield(solinit,'y')
  error(['The field ''y'' not present in ''' inputname(3) ''''])
end

x = solinit.x(:)';     % row vector
y = solinit.y;

if isempty(x) | (length(x)<2)
  error(['''' inputname(3) '.x'' must contain at least the two end points'])
else
  N = length(x);      % number of mesh points
end
xdir = sign(x(end)-x(1));
if any( xdir * diff(x) <= 0 )
  error (['The entries in ''' inputname(3) ...
          '.x'' must strictly increase or decrease']);
end

if isempty(y)
  error(['No initial guess provided in ''' inputname(3) '.y'''])
end
if length(y(1,:)) ~= N
  error(['''' inputname(3) '.y'' not consistent with ''' ...
         inputname(3) '.x'''])
end

n = size(y,1);  % size of the DE system
nN = n*N;

% stats
nODEeval = 0;
nBCeval = 0;

% set the options
if nargin<4
  options = [];
end

% parameters
knownPar = varargin;
unknownPar = isfield(solinit,'parameters');
if unknownPar
  npar = length(solinit.parameters(:));
  ExtraArgs = [{solinit.parameters(:)} knownPar];
else
  npar = 0;
  ExtraArgs = knownPar;
end

% Get options and set the defaults
rtol = bvpgetc(options,'RelTol',1e-3);
if (length(rtol) ~= 1) | (rtol<=0)
  error('RelTol must be a positive scalar.');
end
if rtol < 100*eps
  rtol = 100*eps;
  warning(['RelTol has been increased to ' num2str(rtol) '.']);
end
atol = bvpgetc(options,'AbsTol',1e-6);
if length(atol) == 1
  atol = atol(ones(n,1));
else
  if length(atol) ~= n
    error(sprintf(['Solving %s requires a scalar AbsTol, '...
                   'or a vector AbsTol of length %d'],funstring(ode),n));
  end
  atol = atol(:);
end
if any(atol<=0)
  error('AbsTol must be positive')
end

threshold = atol/rtol;

% analytical Jacobians
Fjac = bvpgetc(options,'FJacobian');
BCjac = bvpgetc(options,'BCJacobian');
Nmax = bvpgetc(options,'Nmax',floor(1000/n));
printstats = strcmp(bvpgetc(options,'Stats','off'),'on');

xyVectorized = strcmp(bvpgetc(options,'Vectorized','off'),'on');
if xyVectorized     % vectorized wrt x, y
  vectorized = 2;   % input to numjac
else
  vectorized = 0;
end

if xyVectorized
  yp = feval(ode,x,y,ExtraArgs{:});
  nODEeval = nODEeval + 1;
else
  yp = zeros(n,N);
  for i=1:N
    yp(:,i) = feval(ode,x(i),y(:,i),ExtraArgs{:});
  end
  nODEeval = nODEeval + N;
end

maxNewtIter = 4;
maxProbes = 4;    % weak line search
needGlobJac = true;

done = false;

% THE MAIN LOOP:
while ~done

  if unknownPar
    Y = [y(:);ExtraArgs{1}];
  else
    Y =  y(:);
  end

  [RHS,Fmid,NF] = colloc_eqns(n,x,Y,yp,ode,bc,npar,xyVectorized,ExtraArgs);
  nODEeval = nODEeval + NF;
  nBCeval = nBCeval + 1;

  for iter=1:maxNewtIter
    if needGlobJac
      % setup and factor the global Jacobian
      [dPHIdy,NF,NBC] = colloc_Jac(n,x,Y,yp,ode,bc,Fjac,BCjac,npar,vectorized,ExtraArgs);
      needGlobJac = false;
      nODEeval = nODEeval + NF;
      nBCeval = nBCeval + NBC;
      % explicit row scaling
      scalMatrix = spdiags(1 ./ max(abs(dPHIdy)')',0,nN+npar,nN+npar);
      dPHIdy = scalMatrix  * dPHIdy;
      % check the Jacobian for singularity
      singJac = (condest(dPHIdy) > 1e10);
      if singJac
        warning('An ill-conditioned Jacobian -- the mesh will be halved');
        break
      end
      [L,U,P] = lu(dPHIdy);
      scalMatrix = P * scalMatrix;
    end

    % find the Newton direction,
    delY = U\(L\( scalMatrix * RHS ));
    distY = norm(delY);

    % weak line search with an affine stopping criterion
    lambda = 1;
    for probe=1:maxProbes
      Ynew = Y - lambda*delY;

      y  = reshape(Ynew(1:nN),n,N);
      if unknownPar
        ExtraArgs{1} = Ynew(nN+1:end);
      end

      if xyVectorized
        yp = feval(ode,x,y,ExtraArgs{:});
        nODEeval = nODEeval + 1;
      else
        yp = zeros(n,N);
        for i=1:N
          yp(:,i) = feval(ode,x(i),y(:,i),ExtraArgs{:});
        end
        nODEeval = nODEeval + N;
      end

      [RHS,Fmid,NF] = colloc_eqns(n,x,Ynew,yp,ode,bc,npar,xyVectorized,ExtraArgs);
      nODEeval = nODEeval + NF;
      nBCeval = nBCeval + 1;

      distYnew = norm(U \ (L \ (scalMatrix * RHS)));

      if (distYnew < 0.9*distY)
        break
      else
        lambda = 0.5*lambda;
      end
    end

    Y  = Ynew;   % y, yp, ExtraArgs, and RHS are consistent with Ynew

    needGlobJac = (distYnew > 0.1*distY);

    if distYnew < 0.1*rtol
      break
    end
  end

  if singJac
    NN = 2*N-1;   % halve the mesh
    if NN > Nmax
      msg = sprintf(...
        [ 'Unable to refine the mesh any further -- '...
          'the Jacobian of the collocation equations is singular']);
      error(msg)
    end
    xx = zeros(1,NN);
    xx(1:2:NN) = x;
    xx(2:2:NN-1) = x(1:N-1) + 0.5*diff(x);

    % a new solution profile
    yy = zeros(n,NN);
    yy(:,1:2:NN) = y;
    temp = struct('x',x,'y',y,'yp',yp,'solver','bvp4cc');
    yy(:,2:2:NN-1) = bvpval(temp,xx(2:2:NN-1));
    yyp = zeros(n,NN);
    yyp(:,1:2:NN) = yp;
    if xyVectorized
      yyp(:,2:2:NN-1) = feval(ode,xx(2:2:NN-1),yy(:,2:2:NN-1),ExtraArgs{:});
      nODEeval = nODEeval + 1;
    else
      for i=2:2:NN-1
        yyp(:,i) = feval(ode,xx(i),yy(:,i),ExtraArgs{:});
      end
      nODEeval = nODEeval + (NN-1)/2;
    end

    N = NN;
    nN = n*N;
    x = xx;
    y = yy;
    yp = yyp;
    needGlobJac = true;

  else

    % residual at the midpoints is related to the error
    % in satisfying the collocation equations
    NewtRes = reshape(RHS(n+npar+1:end),n,N-1);
    temp = ( NewtRes * spdiags(1.5./diff(x(:)),0,N-1,N-1)) ./ ...
             max(abs(Fmid),threshold(:,ones(1,N-1)));
    res_mid = sqrt(dot(temp,temp,1));

    [res,NF] = residual(ode,x,y,yp,res_mid,threshold,xyVectorized,unknownPar,ExtraArgs);
    nODEeval = nODEeval + NF;

    if max(res) < rtol
      done = true;
    else
      % redistribute the mesh
      [N,x,y] = new_profile(n,x,y,yp,res,rtol,Nmax);
      if N>Nmax
        msg = sprintf(...
          [ 'Unable to meet the tolerance without using more than %d '...
            'mesh points. \n The last mesh of %d points and ' ...
            'the solution are available in the output argument. \n ',...
            'The maximum residual is %e, while requested accuracy %e \n'],...
          Nmax,length(x),max(res),rtol);
        warning(msg);
        sol.x = x;
        sol.y = y;
        sol.yp = yp;
        sol.solver = 'bvp4cc';
        if unknownPar
          sol.parameters = ExtraArgs{1};
        end
        return
      end
      nN = n*N;

      if xyVectorized
        yp = feval(ode,x,y,ExtraArgs{:});
        nODEeval = nODEeval + 1;
      else
        yp = zeros(n,N);
        for i=1:N
          yp(:,i) = feval(ode,x(i),y(:,i),ExtraArgs{:});
        end
        nODEeval = nODEeval + N;
      end

      needGlobJac = true;
    end
  end
end    % while

% Output
sol.x = x;
sol.y = y;
sol.yp = yp;
sol.solver = 'bvp4cc';
if unknownPar
  sol.parameters = ExtraArgs{1};
end

% Stats
if printstats
  fprintf('The solution was obtained on a mesh of %g points.\n',N);
  fprintf('The maximum residual is %10.3e. \n',max(res));
  fprintf('There were %g calls to the ODE function. \n',nODEeval);
  fprintf('There were %g calls to the BC function. \n',nBCeval);
end

%--------------------------------------------------------------------------

function [Sx,Spx] = Hermite_cubic (xx, xn,yn,ypn,xnp1,ynp1,ypnp1);
%HERMITE_CUBIC Evaluate the cubic Hermite interpolant and its first
%  derivative

h = xnp1 - xn;
slope = (ynp1-yn)/h;
c = (3*slope - 2*ypn - ypnp1)/h;
d = (ypn+ypnp1-2*slope) / h^2;
t = xx - xn;
Sx = ((d*t+c)*t+ypn)*t+yn;
if nargout>1
  Spx = (3*d*t+2*c)*t+ypn;
end

%---------------------------------------------------------------------------

function [Sx,Spx] = interp_Lobatto_points (hnode,diffx,y,yp);
%INTERP_LOBATTO_POINTS  Evaluate the cubic Hermite interpolant and its first
%  derivative at x+hnode*diffx
N = size(y,2);
diffx = diffx(:);  % column vector
diagscal = spdiags(1./diffx,0,N-1,N-1);
slope = (y(:,2:N) - y(:,1:N-1)) * diagscal;
c = (3*slope - 2*yp(:,1:N-1) - yp(:,2:N)) * diagscal;
d = (yp(:,1:N-1)+yp(:,2:N)-2*slope) * (diagscal.^2);

diagscal = spdiags(hnode*diffx,0,diagscal);
d = d*diagscal;

Sx = ((d+c)*diagscal+yp(:,1:N-1))*diagscal+y(:,1:N-1);
Spx = (3*d+2*c)*diagscal+yp(:,1:N-1);

%---------------------------------------------------------------------------

function [res,nfcn] = residual (Fcn, x, y, yp, mid_res, threshold, xyVectorized, ...
                                unknownPar, ExtraArgs)
%RESIDUAL  Compute L2-norm of the residual using 5-point Lobatto quadrature

% Lobatto quadrature
lob4 = (1+sqrt(3/7))/2;
lob2 = (1-sqrt(3/7))/2;
lobw24 = 49/90;
lobw3 = 32/45;

if xyVectorized
  [n,N] = size(y);
  h = diff(x);   % interval lengths
  thresh = threshold(:,ones(1,N-1));

  % the mid-points
  res = lobw3*mid_res.^2;

  % Lobatto L2 points
  xLob = x(1:N-1)+lob2*h;
  [yLob,ypLob] = interp_Lobatto_points(lob2,h,y,yp);
  fLob = feval(Fcn,xLob,yLob,ExtraArgs{:});
  temp = (ypLob - fLob) ./ max(abs(fLob),thresh);
  resLob = dot(temp,temp,1);
  res = res + lobw24*resLob;

  % Lobatto L4 points
  xLob = x(1:N-1)+lob4*h;
  [yLob,ypLob] = interp_Lobatto_points(lob4,h,y,yp);
  fLob = feval(Fcn,xLob,yLob,ExtraArgs{:});
  temp = (ypLob - fLob) ./ max(abs(fLob),thresh);
  resLob = dot(temp,temp,1);
  res = sqrt( abs(h/2) .* (res + lobw24*resLob));

  nfcn = 2;

else  % not vectorized
  N = length(x);
  res = zeros(1,N-1);
  res_sample = zeros(3,1);
  h = diff(x);   % interval lengths
  for i=1:N-1
    % discrete norm of the residual at Lob2
    x2 = x(i) + lob2*h(i);
    [y2,yp2] = Hermite_cubic(x2,x(i),y(:,i),yp(:,i),x(i+1),y(:,i+1),yp(:,i+1));     f2 = feval(Fcn,x2,y2,ExtraArgs{:});
    res_sample(1) = norm( (yp2-f2) ./ max(abs(f2),threshold) );

    % discrete norm of the residual at Lob3, the midpoint
    res_sample(2) = mid_res(i);

    % discrete norm of the residual at Lob4
    x4 = x(i) + lob4*h(i);
    [y4,yp4] = Hermite_cubic(x4,x(i),y(:,i),yp(:,i),x(i+1),y(:,i+1),yp(:,i+1));     f4 = feval(Fcn,x4,y4,ExtraArgs{:});
    res_sample(3) = norm( (yp4-f4) ./ max(abs(f4),threshold) );

    % the residual on [x(i),x(i+1)]
    temp = abs(h(i)/2) * (lobw24*(res_sample(1)^2 + res_sample(3)^2) +...
                      lobw3*res_sample(2)^2);
    res(i) = sqrt(temp);
  end
  nfcn = 2*(N-1);
end


%---------------------------------------------------------------------------

function [NN,xx,yy] = new_profile(n,x,y,yp,res,rtol,Nmax)
%NEW_PROFILE  Redistribute mesh points and approximate the solution

N = length(x);
i1 = find(res>rtol);
i2 = find(res(i1)>100*rtol);
NNmax = N + length(i1) + length(i2);  % upper bound of the mesh size
n = size(y,1);
xx = zeros(1,NNmax);
yy = zeros(n,NNmax);

last_int = N-1;
h = diff(x);
xx(1) = x(1);
yy(:,1) = y(:,1);
NN = 1;
i=1;
removed = 0;
while i<=last_int
  if res(i) > rtol     % introduce points
    if res(i) > 100*rtol
      Ni = 2;
    else
      Ni = 1;
    end
    hi = h(i) / (Ni+1);
    for j=1:Ni
      xx(NN+j) = xx(NN+j-1) + hi;
      yy(:,NN+j) = Hermite_cubic(xx(NN+j),x(i),y(:,i),yp(:,i),...
                           x(i+1),y(:,i+1),yp(:,i+1));
    end
    NN = NN + Ni;
  else
    if i <= last_int-2
      if max(res(i+1:i+2)) < rtol   % try to remove points
        hnew = (h(i)+h(i+1)+h(i+2))/2;
        C1 = res(i)/(h(i)/hnew)^(7/2);
        C2 = res(i+1)/(h(i+1)/hnew)^(7/2);
        C3 = res(i+2)/(h(i+2)/hnew)^(7/2);
        pred_res = max([C1,C2,C3]);

        if pred_res < 0.5 * rtol   % replace 3 intervals with 2
          removed = removed+1;
          xx(NN+1) = xx(NN) + hnew;
          yy(:,NN+1) = Hermite_cubic(xx(NN+1),x(i+1),y(:,i+1),yp(:,i+1),...
                               x(i+2),y(:,i+2),yp(:,i+2));
          NN = NN+1;
          i = i+2;
        end
      end
    end
  end
  NN = NN+1;
  xx(NN) = x(i+1);   % preserve the next mesh point
  yy(:,NN) = y(:,i+1);
  i = i+1;
end
if NN > Nmax  % too many mesh points requested
  xx = x;
  yy = y;     % return the previous solution
else
  % Trim the output
  xx = xx(1:NN);
  yy = yy(:,1:NN);
end

%---------------------------------------------------------------------------

function [Phi,Fmid,nfcn] = colloc_eqns(n, x, Y, F, Fcn, Gbc, npar, xyVectorized, ExtraArgs)
%COLLOC_EQNS  Evaluate the system of collocation equations Phi(Y).
%   The derivative approximated at the midpoints and returned in Fmid is
%   used to estimate the residual.

N = length(x);
nN = n*N;
Phi = zeros(nN+npar,1);

h = diff(x);
Hdiag = spdiags(h(:),0,N-1,N-1);

xip05 = x(1:N-1) + 0.5*h;

y = reshape(Y(1:nN),n,N);

% Boundary conditions
y1 = y(:,1);       % y(a)
yN = y(:,end);     % y(b)
Phi(1:n+npar) = feval(Gbc,y1,yN,ExtraArgs{:});

% Derivative at the midpoints
if xyVectorized

  yip05 = (y(:,2:N) + y(:,1:N-1))/2 - (F(:,2:N)-F(:,1:N-1))*(Hdiag/8);
  Fmid = feval(Fcn,xip05,yip05,ExtraArgs{:});
  nfcn = 1;
else % not vectorized

  Fmid = zeros(n,N-1);
  yip1 = y(:,1);
  Fip1 = F(:,1);
  for i = 1:N-1
    yi = yip1;
    Fi = Fip1;
    yip1 = y(:,i+1);
    Fip1 = F(:,i+1);
    yip05 = (yi+yip1)/2 - (Fip1-Fi)*(h(i)/8);
    Fmid(:,i) = feval(Fcn,xip05(i),yip05,ExtraArgs{:});
  end
  nfcn = N-1;
end

% the Lobatto IIIa formula
Phi(npar+n+1:npar+nN) = reshape( y(:,2:N) - y(:,1:N-1)...
                                 - (F(:,2:N)+4*Fmid+F(:,1:N-1))*(Hdiag/6),...
                                 nN-n,1);

%---------------------------------------------------------------------------

function [dPHIdy,nfcn,nbc] = colloc_Jac(n, x, Y, F, ode, bc, ...
                                        Fjac, BCjac, npar, vectorized, ...
                                        ExtraArgs);
%colloc_Jac  Form the Jacobian of the system of collocation eqns

nfcn = 0;
nbc = 0;
N = length(x);
nN = n*N;
In = eye(n);

ya = Y(1:n);
yb = Y(nN-n+1:nN);

rows = (n+1:2*n);   % define a region of
cols = (1:n);       % the global Jacobian

if npar == 0    % no unknown parameters

  dPHIdy = spalloc(nN,nN,2*nN*n);  % sparse storage

  % BC
  if isempty(BCjac)   % use numerical approx
    [dGdya,dGdyb,nbcCalls] = BCnumjac(bc,ya,yb,n,npar,ExtraArgs);
  else  % use analytical Jacobian
    [dGdya,dGdyb] = feval(BCjac,ya,yb,ExtraArgs{:});
  end
  dPHIdy(1:n,1:n) = dGdya;  % and postpone storing dGdyb

  % Collocation equations
  if isempty(Fjac)  % use numerical approx

    xip1 = x(1);
    yip1 = Y(cols);
    Fip1 = F(:,1);
    threshval = 1e-6;
    threshold = threshval(ones(n,1));
    FAC = [];
    [Jip1,FAC,dummy,dummy,nFcalls] = ...
     numjac(ode,xip1,yip1,Fip1,threshold,FAC,vectorized,[],[],ExtraArgs{:}); 
    nfcn = nfcn+nFcalls;
    nrmJip1 = norm(Jip1,1);

    for i = 1:N-1
      % the current grid point
      xi = xip1;
      yi = yip1;
      Fi = Fip1;
      Ji = Jip1;
      nrmJi = nrmJip1;
      % the next grid point
      xip1 = x(i+1);
      yip1 = Y(cols+n);
      Fip1 = F(:,i+1);
      [Jip1,FAC,dummy,dummy,nFcalls] = ...
      numjac(ode,xip1,yip1,Fip1,threshold,FAC,vectorized,[],[],ExtraArgs{:});
      nfcn = nfcn+nFcalls;
      nrmJip1 = norm(Jip1,1);

      % the midpoint
      hi = xip1-xi;
      xip05 = (xi+xip1)/2;
      if norm(Jip1 - Ji,1) <= 0.25*(nrmJi + nrmJip1)
          twiceJip05 = Ji + Jip1;
      else
          yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
          Fip05 = feval(ode,xip05,yip05,ExtraArgs{:});
          [Jip05,FAC,dummy,dummy,nFcalls] = ...
     numjac(ode,xip05,yip05,Fip05,threshold,FAC,vectorized,[],[],ExtraArgs{:});
          nfcn = nfcn+nFcalls+1;
          twiceJip05 = 2*Jip05;
      end
      % assembly
      dPHIdy(rows,cols) = -(In+hi/6*(Ji+twiceJip05*(In+hi/4*Ji)));
      cols = cols + n;
      dPHIdy(rows,cols) = In-hi/6*(Jip1+twiceJip05*(In-hi/4*Jip1));
      rows = rows+n;   % next equation
    end
  else   % use analytical Jacobian
    xip1 = x(1);
    yip1 = Y(cols);
    Fip1 = F(:,1);
    Jip1 = feval(Fjac,xip1,yip1,ExtraArgs{:});

    for i = 1:N-1
      % the current grid point
      xi = xip1;
      yi = yip1;
      Fi = Fip1;
      Ji = Jip1;
      % the next grid point
      xip1 = x(i+1);
      yip1 = Y(cols+n);
      Fip1 = F(:,i+1);
      Jip1 = feval(Fjac,xip1,yip1,ExtraArgs{:});
      % the midpoint
      hi = xip1-xi;
      xip05 = (xi+xip1)/2;
      yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
      Jip05 = feval(Fjac,xip05,yip05,ExtraArgs{:});  % recompute the Jacobian
      % assembly
      twiceJip05 = 2*Jip05;
      dPHIdy(rows,cols) = -(In+hi/6*(Ji+twiceJip05*(In+hi/4*Ji)));
      cols = cols + n;
      dPHIdy(rows,cols) = In-hi/6*(Jip1+twiceJip05*(In-hi/4*Jip1));
      rows = rows+n;   % next equation
    end
  end
  dPHIdy(1:n,nN-n+1:nN) = dGdyb;    % computed earlier

else  % there are unknown parameters

  dPHIdy = spalloc(nN+npar,nN+npar,(nN+npar)*(2*n+npar));  % sparse storage


  last_cols = zeros(nN+npar,npar);   % accumulator

  % BC
  if isempty(BCjac)   % use numerical approx
    [dGdya,dGdyb,nbcCalls,dGdpar] = BCnumjac(bc,ya,yb,n,npar,ExtraArgs);
  else  % use analytical Jacobians
    [dGdya,dGdyb,dGdpar] = feval(BCjac,ya,yb,ExtraArgs{:});
  end
  dPHIdy(1:n+npar,1:n) = dGdya;  % and postpone storing dGdyb
  last_cols(1:n+npar,:) = dGdpar;

  % accomodate additional npar of boundary conditions
  rows = rows+npar;

  if isempty(Fjac)  % use numerical approx
    xip1 = x(1);
    yip1 = Y(cols);
    Fip1 = F(:,1);
    threshval = 1e-6;
    threshold = threshval(ones(n,1));
    FAC = [];
    threshpar = threshval(ones(npar,1));
    FACpar = [];
    [Jip1,FAC,dFdpar_ip1,FACpar,nFcalls] = ...
        Fnumjac(ode,xip1,yip1,Fip1,threshold,FAC,vectorized,...
                                   threshpar,FACpar,ExtraArgs);
    nfcn = nfcn+nFcalls;
    nrmJip1 = norm(Jip1,1);
    nrmdFdpar_ip1 = norm(dFdpar_ip1,1);
    for i = 1:N-1
      % the current grid point
      xi = xip1;
      yi = yip1;
      Fi = Fip1;
      Ji = Jip1;
      dFdpar_i = dFdpar_ip1;
      nrmJi = nrmJip1;
      nrmdFdpar_i = nrmdFdpar_ip1;
      % the next grid point
      xip1 = x(i+1);
      yip1 = Y(cols+n);
      Fip1 = F(:,i+1);
      [Jip1,FAC,dFdpar_ip1,FACpar,nFcalls] = ...
          Fnumjac(ode,xip1,yip1,Fip1,threshold,FAC,vectorized,...
                                     threshpar,FACpar,ExtraArgs);
      nfcn = nfcn+nFcalls;
      nrmJip1 = norm(Jip1,1);
      nrmdFdpar_ip1 = norm(dFdpar_ip1,1);
      % the midpoint
      hi = xip1-xi;
      xip05 = (xi+xip1)/2;
      if (norm(Jip1 - Ji,1) <= 0.25*(nrmJi + nrmJip1)) & ...
         (norm(dFdpar_ip1 - dFdpar_i,1) <= 0.25*(nrmdFdpar_i + nrmdFdpar_ip1))
        Jip05 = 0.5*(Ji + Jip1);
        dFdpar_ip05 = 0.5*(dFdpar_i + dFdpar_ip1);
      else
        yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
        Fip05 = feval(ode,xip05,yip05,ExtraArgs{:});
        [Jip05,FAC,dFdpar_ip05,FACpar,nFcalls] = ...
            Fnumjac(ode,xip05,yip05,Fip05,threshold,FAC,vectorized,...
                                         threshpar,FACpar,ExtraArgs);

        nfcn = nfcn+nFcalls+1;
      end
      twiceJip05 = 2*Jip05;
      % assembly
      dPHIdy(rows,cols) = -(In+hi/6*(Ji+twiceJip05*(In+hi/4*Ji)));
      cols = cols + n;
      dPHIdy(rows,cols) = In-hi/6*(Jip1+twiceJip05*(In-hi/4*Jip1));
      last_cols(rows,:) = -hi*dFdpar_ip05 + hi^2/12*Jip05*...
                                            (dFdpar_ip1-dFdpar_i);
      rows = rows+n;   % next equation
    end
  else   % use analytical Jacobians
    xip1 = x(1);
    yip1 = Y(cols);
    Fip1 = F(:,1);
    [Jip1,dFdpar_ip1] = feval(Fjac,xip1,yip1,ExtraArgs{:});
    for i = 1:N-1
      % the current grid point
      xi = xip1;
      yi = yip1;
      Fi = Fip1;
      Ji = Jip1;
      dFdpar_i = dFdpar_ip1;
      % the next grid point
      xip1 = x(i+1);
      yip1 = Y(cols+n);
      Fip1 = F(:,i+1);
      [Jip1, dFdpar_ip1] = feval(Fjac,xip1,yip1,ExtraArgs{:});
      % the midpoint
      hi = xip1-xi;
      xip05 = (xi+xip1)/2;
      yip05 = (yi+yip1)/2-hi/8*(Fip1-Fi);
      [Jip05, dFdpar_ip05] = feval(Fjac,xip05,yip05,ExtraArgs{:});
      % assembly
      twiceJip05 = 2*Jip05;
      dPHIdy(rows,cols) = -(In+hi/6*(Ji+twiceJip05*(In+hi/4*Ji)));
      cols = cols+n;
      dPHIdy(rows,cols) = In-hi/6*(Jip1+twiceJip05*(In-hi/4*Jip1));
      last_cols(rows,:) = -hi*dFdpar_ip05 + hi^2/12*Jip05* ...
                                            (dFdpar_ip1-dFdpar_i);
      rows = rows+n;   % next equation
    end
  end
  dPHIdy(1:n+npar,nN-n+1:nN) = dGdyb;    % computed earlier
  dPHIdy(:,nN+1:nN+npar) = last_cols;  % accumulated

end

ndfdy = 2*N-1;

%---------------------------------------------------------------------------

function [dBCdya,dBCdyb,nbcCalls,dBCdpar] = BCnumjac(bc,ya,yb,n,npar,ExtraArgs)
if npar > 0   % unknown parameters
  argvect = [ya; yb; ExtraArgs{1}];   % stack column vectors
else
  argvect = [ya; yb];        % stack column vectors
end
BCval = feval(bc,ya,yb,ExtraArgs{:});
threshval = 1e-6;
thresh = threshval(ones(2*n+npar,1));
fac = [];
vectorized = 0; % do not vectorize BC function
[temp,fac,dummy,dummy,nbcCalls] = numjac(@BCaux,bc,argvect,BCval,...
                                  thresh,fac,vectorized,[],[],n,npar,ExtraArgs);dBCdya = temp(:,1:n);
dBCdyb = temp(:,n+1:2*n);
if npar>0  % unknown parameters
  dBCdpar = temp(:,2*n+1:end);
end
nbcCalls = nbcCalls+1; % stats

%---------------------------------------------------------------------------

function BCval = BCaux(bc,argvec,n,npar,ExtraArgs);
ya = argvec(1:n);
yb = argvec(n+1:2*n);
if npar > 0
  ExtraArgs{1} = argvec(2*n+1:end);  % unknown parameters
end
BCval = feval(bc,ya,yb,ExtraArgs{:});
 
%---------------------------------------------------------------------------
 
function [Ji,FAC,dFdpari,FACpar,nFcalls] = Fnumjac(ode,xi,yi,Fi,thresh,FAC,vectorized,...
threshpar,FACpar,ExtraArgs);
%Compute dF/dy and dF/dpar
[Ji,FAC,dummy,dummy,nFdyCalls] = numjac(ode,xi,yi,Fi,thresh,FAC,...
                                  vectorized,[],[],ExtraArgs{:});
pari = ExtraArgs{1};   % unknown parameters
[dFdpari,FACpar,dummy,dummy,nFdparCalls] = numjac(@Faux,xi,pari,Fi,threshpar,FACpar,...
                                             0,[],[],ode,yi,ExtraArgs);
nFcalls = nFdyCalls + nFdparCalls;
 
%---------------------------------------------------------------------------
 
function Fval = Faux(x,par,ode,y,ExtraArgs);
ExtraArgs{1} = par;    % unknown parameters
Fval = feval(ode,x,y,ExtraArgs{:});

function o = bvpgetc(options,name,default)
%BVPGET  Get BVP OPTIONS parameters.
%   VAL = BVPGET(OPTIONS,'NAME') extracts the value of the named property
%   from integrator options structure OPTIONS, returning an empty matrix if
%   the property value is not specified in OPTIONS. It is sufficient to type
%   only the leading characters that uniquely identify the property. Case is
%   ignored for property names. [] is a valid OPTIONS argument.
%
%   VAL = BVPGET(OPTIONS,'NAME',DEFAULT) extracts the named property as
%   above, but returns VAL = DEFAULT if the named property is not specified
%   in OPTIONS. For example
%
%       val = bvpget(opts,'RelTol',1e-4);
%
%   returns val = 1e-4 if the RelTol property is not specified in opts.
%
%   See also BVPSET, BVPINIT, BVP4C, DEVAL.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2001 The MathWorks, Inc.
%   $Revision: 1.8 $  $Date: 2001/02/16 16:19:39 $

if nargin < 2
  error('Not enough input arguments.');
end
if nargin < 3
  default = [];
end

if ~isempty(options) & ~isa(options,'struct')
  error('First argument must be an options structure created with BVPSET.');
end
 
if isempty(options)
  o = default;
  return;
end
 
Names = [
    'AbsTol    '
    'RelTol    '
    'FJacobian '
    'BCJacobian'
    'Stats     '
    'Nmax      '
    'Vectorized'
    ];
 
[m,n] = size(Names);
names = lower(Names);
 
lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j)               % if no matches
  error(sprintf(['Unrecognized property name ''%s''.  ' ...
                 'See BVPSET for possibilities.'], name));
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    msg = sprintf('Ambiguous property name ''%s'' ', name);
    msg = [msg '(' deblank(Names(j(1),:))];
    for k = j(2:length(j))'
      msg = [msg ', ' deblank(Names(k,:))];
    end
    msg = sprintf('%s).', msg);
    error(msg);
  end
end
 
if any(strcmp(fieldnames(options),deblank(Names(j,:))))
  o = getfield(options, Names(j,:));
  if isempty(o)
    o = default;
  end
else
  o = default;
end

function options = bvpsetc(varargin)
%BVPSET  Create/alter BVP OPTIONS structure.
%   OPTIONS = BVPSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. It is
%   sufficient to type only the leading characters that uniquely identify the
%   property. Case is ignored for property names.
%
%   OPTIONS = BVPSET(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS.
%
%   OPTIONS = BVPSET(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties overwrite
%   corresponding old properties.
%
%   BVPSET with no input arguments displays all property names and their
%   possible values.
%
%BVPSET PROPERTIES
%
%RelTol - Relative tolerance for the residual [ positive scalar {1e-3} ]
%   This scalar applies to all components of the residual vector, and
%   defaults to 1e-3 (0.1% accuracy). The computed solution S(x) is the exact
%   solution of S'(x) = F(x,S(x)) + res(x). On each subinterval of the mesh,
%   component i of the residual satisfies
%          norm( res(i) / max( [abs(F(i)) , AbsTol(i)/RelTol] ) ) <= RelTol.
%
%AbsTol - Absolute tolerance for the residual [ positive scalar or vector {1e-6} ]
%   A scalar tolerance applies to all components of the residual vector.
%   Elements of a vector of tolerances apply to corresponding components of
%   the residual vector. AbsTol defaults to 1e-6. See RelTol.
%
%FJacobian - Analytical partial derivatives of ODEFUN [ function ]
%   For example, when solving y' = f(x,y), set this property to @FJAC if
%   DFDY = FJAC(X,Y) evaluates the Jacobian of f with respect to y.
%   If the problem involves unknown parameters, [DFDY,DFDP] = FJAC(X,Y,P)
%   must also return the partial derivative of f with respect to p.
%
%BCJacobian - Analytical partial derivatives of BCFUN [ function ]
%   For example, for boundary conditions bc(ya,yb) = 0, set this property to
%   @BCJAC if [DBCDYA,DBCDYB] = BCJAC(YA,YB) evaluates the partial
%   derivatives of bc with respect to ya and to yb. If the problem involves
%   unknown parameters, [DBCDYA,DBCDYB,DBCDP] = BCJAC(YA,YB,P) must also
%   return the partial derivative of bc with respect to p.
%
%Nmax - Maximum number of mesh points allowed [positive integer {floor(1000/n)}]%
%Stats - Display computational cost statistics  [ on | {off} ]
%
%Vectorized - Vectorized ODE function  [ on | {off} ]
%   Set this property 'on' if the derivative function
%   ODEFUN([x1 x2 ...],[y1 y2 ...]) returns [ODEFUN(x1,y1) ODEFUN(x2,y2) ...].
%   When parameters are present, the derivative function
%   ODEFUN([x1 x2 ...],[y1 y2 ...],p) should return
%   [ODEFUN(x1,y1,p) ODEFUN(x2,y2,p) ...].
%
%   See also BVPGET, BVPINIT, BVP4C, DEVAL.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2001 The MathWorks, Inc.
%   $Revision: 1.7 $  $Date: 2001/02/16 16:19:43 $

% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('          AbsTol: [ positive scalar or vector {1e-6} ]\n');
  fprintf('          RelTol: [ positive scalar {1e-3} ]\n');
  fprintf('       FJacobian: [ function ]\n');
  fprintf('      BCJacobian: [ function ]\n');
  fprintf('           Stats: [ on | {off} ]\n');
  fprintf('            Nmax: [ nonnegative integer {floor(1000/n)} ]\n');
  fprintf('      Vectorized: [ on | {off} ]\n');
  fprintf('\n');
  return;
end

Names = [
    'AbsTol    '
    'RelTol    '
    'FJacobian '
    'BCJacobian'
    'Stats     '
    'Nmax      '
    'Vectorized'
    ];

[m,n] = size(Names);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
options = [];
i = 1;
while i <= nargin
  arg = varargin{i};
  if isstr(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(sprintf(['Expected argument %d to be a string property name '...
                     'or an options structure\ncreated with BVPSET.'], i));
    end
    if isempty(options)
      options = arg;
    else
      for j = 1:m
        val = getfield(arg,Names(j,:));
        if ~isequal(val,[])             % empty strings '' do overwrite
          options = setfield(options,Names(j,:),val);
        end
      end
    end
  end
  i = i + 1;
end
if isempty(options)
  for j = 1:m
    options = setfield(options,Names(j,:),[]);
  end
end
 
% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end
expectval = 0;                      % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
  if ~expectval
    if ~isstr(arg)
      error(...
        sprintf('Expected argument %d to be a string property name.', i));
    end
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(sprintf('Unrecognized property name ''%s''.', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error(msg);
      end
    end
    expectval = 1;                      % we expect a value next
  else
    options = setfield(options,Names(j,:),arg);
    expectval = 0;
  end
  i = i + 1;
end
 
if expectval
  error(sprintf('Expected value for property ''%s''.', arg));
end
function [dFdy,fac,g,nfevals,nfcalls] = ...
    numjac(F,t,y,Fty,thresh,fac,vectorized,S,g,varargin)
%NUMJAC Numerically compute the Jacobian dF/dY of function F(T,Y).
%   [DFDY,FAC] = NUMJAC(F,T,Y,FTY,THRESH,FAC,VECTORIZED) numerically
%   computes the Jacobian of function F(T,Y), returning it as full matrix
%   DFDY.  T is the independent variable and column vector Y contains the
%   dependent variables. Function F must return a column vector. Vector FTY
%   is F evaluated at (T,Y). Column vector THRESH provides a threshold of
%   significance for Y, i.e. the exact value of a component Y(i) with
%   abs(Y(i)) < THRESH(i) is not important. All components of THRESH must
%   be positive. Column FAC is working storage. On the first call, set
%   FAC to []. Do not alter the returned value between calls. VECTORIZED
%   tells NUMJAC whether multiple values of F can be obtained with a single
%   function evaluation. In particular, VECTORIZED=1 indicates that
%   F(t,[y1 y2 ...]) returns [F(t,y1) F(t,y2) ...], and VECTORIZED=2 that
%   F([x1 x2 ...],[y1 y2 ...]) returns [F(x1,y1) F(x2,y2) ...]. When solving
%   ODEs, use ODESET to set the ODE solver 'Vectorized' property to 'on' if
%   the ODE function has been coded so that F(t,[y1 y2 ...]) returns
%   [F(t,y1) F(t,y2) ...]. When solving BVPs, use BVPSET to set the BVP
%   solver 'Vectorized' property to 'on' if the ODE function has been coded
%   so that F([x1 x2 ...],[y1 y2 ...]) returns [F(x1,y1) F(x2,y2) ...].
%   Vectorizing the function F may speed up the computation of DFDY.
%
%   [DFDY,FAC,G] = NUMJAC(F,T,Y,FTY,THRESH,FAC,VECTORIZED,S,G) numerically
%   computes a sparse Jacobian matrix DFDY.  S is a non-empty sparse matrix
%   of 0's and 1's.  A value of 0 for S(i,j) means that component i of the
%   function F(T,Y) does not depend on component j of vector Y (hence
%   DFDY(i,j) = 0).  Column vector G is working storage.  On the first call,
%   set G to [].  Do not alter the returned value between calls.
%
%   [DFDY,FAC,G,NFEVALS,NFCALLS] = NUMJAC(...) returns the number of values
%   of computed while forming dFdy (NFEVALS) and the number of calls to the
%   function F (NFCALLS). If F is not vectorized, NFCALLS equals NFEVALS.
%
%   Although NUMJAC was developed specifically for the approximation of
%   partial derivatives when integrating a system of ODE's, it can be used
%   for other applications.  In particular, when the length of the vector
%   returned by F(T,Y) is different from the length of Y, DFDY is
%   rectangular.
%
%   See also COLGROUP, ODE15S, ODE23S, ODE23T, ODE23TB, ODESET.

%   NUMJAC is an implementation of an exceptionally robust scheme due to
%   Salane for the approximation of partial derivatives when integrating
%   a system of ODEs, Y' = F(T,Y). It is called when the ODE code has an
%   approximation Y at time T and is about to step to T+H.  The ODE code
%   controls the error in Y to be less than the absolute error tolerance
%   ATOL = THRESH.  Experience computing partial derivatives at previous
%   steps is recorded in FAC.  A sparse Jacobian is computed efficiently
%   by using COLGROUP(S) to find groups of columns of DFDY that can be
%   approximated with a single call to function F.  COLGROUP tries two
%   schemes (first-fit and first-fit after reverse COLMMD ordering) and
%   returns the better grouping.
%
%   D.E. Salane, "Adaptive Routines for Forming Jacobians Numerically",
%   SAND86-1319, Sandia National Laboratories, 1986.
%
%   T.F. Coleman, B.S. Garbow, and J.J. More, Software for estimating
%   sparse Jacobian matrices, ACM Trans. Math. Software, 11(1984)
%   329-345.
%
%   L.F. Shampine and M.W. Reichelt, The MATLAB ODE Suite, SIAM Journal on
%   Scientific Computing, 18-1, 1997.

%   Mark W. Reichelt and Lawrence F. Shampine, 3-28-94
%   Copyright 1984-2001 The MathWorks, Inc.
%   $Revision: 1.34 $  $Date: 2001/04/15 11:59:21 $

% Deal with missing arguments.
if nargin < 10
  args = {};                          % F accepts only (t,y)

  if nargin == 7
    S = [];
  elseif nargin == 6
    S = [];
    vectorized = 0;
  elseif nargin == 5
    S = [];
    vectorized = 0;
    fac = [];
  end
else
  args = varargin;
end

% Initialize.
br = eps ^ (0.875);
bl = eps ^ (0.75);
bu = eps ^ (0.25);
facmin = eps ^ (0.78);
facmax = 0.1;
ny = length(y);
nF = length(Fty);
if isempty(fac)
  fac = sqrt(eps) + zeros(ny,1);
end

% Select an increment del for a difference approximation to
% column j of dFdy.  The vector fac accounts for experience
% gained in previous calls to numjac.
yscale = max(abs(y),thresh);
del = (y + fac .* yscale) - y;
for j = find(del == 0)'
  while 1
    if fac(j) < facmax
      fac(j) = min(100*fac(j),facmax);
      del(j) = (y(j) + fac(j)*yscale(j)) - y(j);
      if del(j)
        break
      end
    else
      del(j) = thresh(j);
      break;
    end
  end
end
if nF == ny
  s = (sign(Fty) >= 0);
  del = (s - (~s)) .* abs(del);         % keep del pointing into region
end

% Form a difference approximation to all columns of dFdy.
if isempty(S)                           % generate full matrix dFdy
  g = [];
  ydel = y(:,ones(1,ny)) + diag(del);
  switch vectorized
   case 1
    Fdel = feval(F,t,ydel,args{:});
    nfcalls = 1;                        % stats
   case 2
    Fdel = feval(F,t(ones(1,ny)),ydel,args{:});
    nfcalls = 1;                        % stats
   otherwise % not vectorized
    Fdel = zeros(nF,ny);
    for j = 1:ny
      Fdel(:,j) = feval(F,t,ydel(:,j),args{:});
    end
    nfcalls = ny;                       % stats
  end
  nfevals = ny;                         % stats (at least one per loop)
  Fdiff = Fdel - Fty(:,ones(1,ny));
  dFdy = Fdiff * diag(1 ./ del);
  [Difmax,Rowmax] = max(abs(Fdiff),[],1);
  % If Fdel is a column vector, then index is a scalar, so indexing is okay.
  absFdelRm = abs(Fdel((0:ny-1)*nF + Rowmax));
else                                    % sparse dFdy with structure S
  if isempty(g)
    g = colgroup(S);                    % Determine the column grouping.
  end
  ng = max(g);
  one2ny = (1:ny)';
  ydel = y(:,ones(1,ng));
  i = (g-1)*ny + one2ny;
  ydel(i) = ydel(i) + del;
  switch vectorized
   case 1
    Fdel = feval(F,t,ydel,args{:});
    nfcalls = 1;                        % stats
   case 2
    Fdel = feval(F,t(ones(1,ng)),ydel,args{:});
    nfcalls = 1;                        % stats
   otherwise % not vectorized
    Fdel = zeros(nF,ng);
    for j = 1:ng
      Fdel(:,j) = feval(F,t,ydel(:,j),args{:});
    end
    nfcalls = ng;                       % stats
  end
  nfevals = ng;                         % stats (at least one per column)
  Fdiff = Fdel - Fty(:,ones(1,ng));
  [i j] = find(S);
  Fdiff = sparse(i,j,Fdiff((g(j)-1)*nF + i),nF,ny);
  dFdy = Fdiff * sparse(one2ny,one2ny,1 ./ del,ny,ny);
  [Difmax,Rowmax] = max(abs(Fdiff),[],1);
  Difmax = full(Difmax);
  % If ng==1, then Fdel is a column vector although index may be a row vector.
  absFdelRm = abs(Fdel((g-1)*nF + Rowmax').');
end

% Adjust fac for next call to numjac.
absFty = abs(Fty);
absFtyRm = absFty(Rowmax);              % not a col vec if absFty scalar
absFtyRm = absFtyRm(:)';                % ensure that absFtyRm is a row vector
j = ((absFdelRm ~= 0) & (absFtyRm ~= 0)) | (Difmax == 0);
if any(j)
  ydel = y;
  Fscale = max(absFdelRm,absFtyRm);

  % If the difference in f values is so small that the column might be just
  % roundoff error, try a bigger increment.
  k1 = (Difmax <= br*Fscale);           % Difmax and Fscale might be zero
  for k = find(j & k1)
    tmpfac = min(sqrt(fac(k)),facmax);
    del = (y(k) + tmpfac*yscale(k)) - y(k);
    if (tmpfac ~= fac(k)) & (del ~= 0)
      if nF == ny
        if Fty(k) >= 0                  % keep del pointing into region
          del = abs(del);
        else
          del = -abs(del);
        end
      end

      ydel(k) = y(k) + del;
      fdel = feval(F,t,ydel,args{:});
      nfevals = nfevals + 1;            % stats
      nfcalls = nfcalls + 1;            % stats
      ydel(k) = y(k);
      fdiff = fdel - Fty;
      tmp = fdiff ./ del;
 
      [difmax,rowmax] = max(abs(fdiff));
      if tmpfac * norm(tmp,inf) >= norm(dFdy(:,k),inf);
        % The new difference is more significant, so
        % use the column computed with this increment.
        if isempty(S)
          dFdy(:,k) = tmp;
        else
          i = find(S(:,k));
          if ~isempty(i)
            dFdy(i,k) = tmp(i);
          end
        end
 
        % Adjust fac for the next call to numjac.
        fscale = max(abs(fdel(rowmax)),absFty(rowmax));
 
        if difmax <= bl*fscale
          % The difference is small, so increase the increment.
          fac(k) = min(10*tmpfac, facmax);
 
        elseif difmax > bu*fscale
          % The difference is large, so reduce the increment.
          fac(k) = max(0.1*tmpfac, facmin);
 
        else
          fac(k) = tmpfac;
 
        end
      end
    end
  end
 
  % If the difference is small, increase the increment.
  k = find(j & ~k1 & (Difmax <= bl*Fscale));
  if ~isempty(k)
    fac(k) = min(10*fac(k), facmax);
  end
 
  % If the difference is large, reduce the increment.
  k = find(j & (Difmax > bu*Fscale));
  if ~isempty(k)
    fac(k) = max(0.1*fac(k), facmin);
  end
end

 
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

res = [  ya(1)-1+1e-4
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

