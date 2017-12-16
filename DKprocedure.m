% This function performs multiview triangulation and translation estimation
% in the case of known rotations. The matrices of internal parameters are
% also assumed to be known.

% Input
%    M - a 4xn matrix, where n is the number of 2D image points. 
%        Each column of M contains the (x,y) coordiantes of the image
%        point, the index of the 3D point of which it is the image and the
%        index of the camera that has recorded this image point. Example:
%        if the jth column of M is (5,156,254,3)', it means that the scene
%        point number 254 has been recorded by the camera number 3 and the
%        resulting image point has coordiantes (5,156).
%    P - a cell of Ncam objects, where Ncam is the total number of cameras.
%        Each element of this cell is either a 3x3 or a 3x4 matrix. In both
%        cases, only the first three columns of this matrices are used.
%        These columns encapsulate the rotations and the internal
%        calibrations.
%    sigma - a positive real number corresponding to the presumed maximal
%        measurement error for inliers. In most cases, sigma=0.5 pixel
%        works well.
%   eps - a positive number used in the bisection algorithm. It should be
%        much smaller than sigma. The standard choice is eps=0.1. A smaller
%        eps will result in larger execution time and in a more accurate
%        approximation of the L_inf cost minimizer. It however does not
%        mean that the resulting estimator of scene points and camera 
%        centers will be closer to the true ones.

% Output
%    Scenepoints - a matrix each column of which contains the coordinates 
%                  of the estimated 3D scene points.
%    CamPos      - a matrix each column of which contains the coordinates 
%                  of the estimated camera centers.
%    Outliers    - a logical vector of the size equal to the number of 2D
%                  measurements. If Outliers(i)=1, then the ith measurement
%                  is an outlier (according to our algorithm).
%    finalRE     - a positive number representing the maximal reprojection 
%                  error of the obtained estimator.

function [Scenepoints,CamPos,Outliers,finalRE]=DKprocedure(M,P,sigma,eps);

    t = clock; 
    [A,C]=ComputeAC(M,P);

    [theta,indexOut]=MyWeightedEstimator(M,A,C,sigma,eps);

    Ncam=max(M(4,:));
    N=2*length(indexOut);

    m(1:2:N)=indexOut;
    m((1:2:N)+1)=indexOut;
    m11=sum(A(:,~m)'~=0)>=4;
    m11(1:3:end)=(m11(1:3:end)|m11(2:3:end)|m11(3:3:end));
    m11(2:3:end)=m11(1:3:end);
    m11(3:3:end)=m11(1:3:end);

    Nout3D=sum(~m11)/3;
    disp([' '])    
    disp(['  ***************************'])
    disp(['  *   removed 2D points: ' num2str(sum(indexOut))]);
    disp(['  *   removed 3D points: ' num2str(Nout3D)]);


    
    theta=theta(m11);
    l=length(theta);
    ind1=l-(3*(Ncam-1));
    Scenepoints=reshape(theta(1:ind1),3,(l/3)-Ncam+1);
    CamPos=(-1)*reshape(theta((l-3*(Ncam-1)+1):l),3,Ncam-1);
    CamPos=[zeros(3,1) CamPos];
    Outliers=indexOut;

    A=A(m11,~m);
    C=C(m11,~m);
    
    finalRE=max(abs(A'*theta)./(C'*theta));
    time=etime(clock,t);
    disp(['  *   Elapsed time is ' num2str(time) ' seconds']);
    disp(['  ***************************'])


%**************************************************************************
%
%     This function takes as input the measurement matrix M1, the matrices
%     A and C, as well as the measurement error sigma and the precision
%     parameter of the bisection algorithm, epsilon. The output is a 
%     This is the last version 
%     July 28, 2009
%
%**************************************************************************



function [theta,indexOut]=MyWeightedEstimator(M1,A,C,ssigma,eps);

[M,N]=size(A);
[womega,wtheta]=WeightedOutlierRemoval(A,C,ssigma);
m=sparse((abs(womega)./(C'*wtheta)>=ssigma/4));
mm=(m(1:2:end)+m(2:2:end))>0; 
clear m;



mmm=((1:N)>0);
mmm(1:2:end)=mm;
mmm(2:2:end)=mm;

%clear mm
A=A(:,~mmm);
C=C(:,~mmm);

clear mmm
[theta]=BisectionLP(A,C,eps,2*max(abs(A'*wtheta)./(C'*wtheta)),0);

indexOut=mm;



%**************************************************************************
%
%     This function forms the matrices A and C from the measurement matrix
%     M1 and the camera matrix P1. Here, P1 is a cell each element of which
%     is a 3x3 or 3x4 camera matrix. In the case when the elements of P1
%     are of size 3x4, only the first three columns are used.
%     This is the last version 
%     July 28, 2009
%
%**************************************************************************



function [A,C]=ComputeAC(M1,P1);

Ncam=max(M1(4,:));
Nobs=2*size(M1,2);
N3D=max(M1(3,:));

A=SPALLOC(3*(Ncam+N3D-1),Nobs,6*Nobs);
C=A;

for (p=1:(Nobs/2))
    x=M1(1:2,p);
    Point3D=M1(3,p);
    Cam=M1(4,p);
    if (Cam==0)
        Cam=Ncam;
    end;
    A((3*Point3D-2):(3*Point3D),(2*p-1):(2*p))=(x*P1{Cam}(3,1:3)-P1{Cam}(1:2,1:3))';
    C((3*Point3D-2):(3*Point3D),(2*p-1):(2*p))=([1;1]*P1{Cam}(3,1:3))';    
    if (Cam>1)
        A((3*N3D+3*(Cam-2)+1):(3*N3D+3*(Cam-1)),(2*p-1):(2*p))=(x*P1{Cam}(3,1:3)-P1{Cam}(1:2,1:3))';
        C((3*N3D+3*(Cam-2)+1):(3*N3D+3*(Cam-1)),(2*p-1):(2*p))=([1;1]*P1{Cam}(3,1:3))';
    end
end;
A=sparse(A);
C=sparse(C);


%**************************************************************************
%
%     Outlier removal by linear programming: the case of weighted L_1-norm
%     This is the last version 
%     July 28, 2009
%
%**************************************************************************

function [womega,wtheta]=WeightedOutlierRemoval(A,C,ssigma);

[M,N]=size(A);
beta=zeros(2*N+M,1);
b=SPALLOC(1,2*N+M,N);
b=sparse((-1)*[ones(1,N) zeros(1,N+M)]');
I=SPEYE(N,N);

% The constraints of the type (omega_j <= s_j) and (-omega_j <= s_j)  
Constraint1=SPALLOC(2*N,2*N+M,4*N); 
Constraint1=[(-1)*I I SPALLOC(N,M,0);(-1)*I (-1)*I SPALLOC(N,M,0)];

% The constraints of the type (|a_p^T*theta-omega_p|<= ssigma*c_p^T*theta)
Constraint21=SPALLOC(2*N,2*N+M,7*N);
Constraint21=SPALLOC(N,2*N+M,7*N);
Constraint21=[SPALLOC(N,N,0) I (-1)*(A'+ssigma*C')];
Constraint22=[SPALLOC(N,N,0) (-1)*I (A'-ssigma*C')];
% The constraint of identifiability: the minimal depth  
% is set to one.

Constraint3=[SPALLOC(N,2*N,0) (-1)*C'];
% AT is the transpose of the matrix appearing in linear constraints.
AT=[Constraint1;Constraint21;Constraint22;Constraint3];

clear I Constraint1 Constraint21 Constraint22 Constraint3; 

% If beta is the variable of minimization, the constraints have the
% form AT*beta <= cc (inequality is understood componentwise).
cc=sparse([SPALLOC(4*N,1,0); (-1)*ones(N,1)]);


disp(['   ------------------------------']);
disp(['   The program starts the fisrt optimization']);


pars.eps=10^(-4);

K.l=length(cc);
pars.fid=0;
[x,beta,info]=sedumi(AT',b,cc,K,pars);


% omega is the vector of outliers: if w_j = 0, then the jth observation
% is an inlier.
omega=beta((N+1):(2*N));

% theta contains the estimator of the pose and the structure
theta=beta((2*N+1):(2*N+M));

% 
%  To enhance the accuracy of estimation, we perform a weighted
%  minimization with wights proportional to the inverse of the depth.
%

disp(['   ------------------------------']);
disp(['   The program starts the weighted optimization']);


b=sparse((-1)*[1./(C'*theta); zeros(N+M,1)]);
pars.eps=10^(-4);
K.l=length(cc);

[x,beta,info]=sedumi(AT',b,cc,K,pars);
disp(['   ------------------------------']);
disp(['   The program starts the bisection algorithm']);
disp(['   ------------------------------']);
womega=beta((N+1):(2*N));
wtheta=beta((2*N+1):(2*N+M));





%**************************************************************************
%
% Bisection algorithm for L_\infty cost minimization under chairality
% constraints
%
%**************************************************************************


function [tttheta]=BisectionLP(A,C,eps,gammaUp,gammaLow);

[M,N]=size(A);
k=0;
while (gammaUp-gammaLow>eps)
    k=k+1;
    gamma=(gammaUp+gammaLow)/2;
    % The constraints of the type (|a_p^T*theta|<= gamma*c_p^T*theta)
    Constraint21=[(-1)*(A'+gamma*C')];
    Constraint22=[A'-gamma*C'];

    % The constraint of chairality
    Constraint3=[(-1)*C'];
    
    % AT is the transpose of the matrix appearing in linear constraints.
    AT=[Constraint21;Constraint22;Constraint3];

    clear Constraint21 Constraint22 Constraint3 ; 

    % If beta is the variable of minimization, the constraints have the
    % form AT*beta <= cc (inequality is understood componentwise).
    cc=sparse([zeros(2*N,1);(-1)*ones(N,1)]);
    pars.eps=10^(-4);
    K.l=length(cc);
    pars.fid=0;
    [x,theta1,info]=sedumi(AT',0,cc,K,pars);
 
    if (info.dinf==1)
        gammaLow=gamma;
    else
        gammaUp1=max(abs(A'*theta1)./(C'*theta1));
        if gammaUp1<=gamma
            tttheta=theta1;
            gammaUp=gammaUp1;
        else
            gammaLow=gamma;
        end;
    end;
    disp(['                   iter ' num2str(k) ': gamma is in the interval [' num2str(gammaLow) ', ' num2str(gammaUp) ']' ]);
end;
