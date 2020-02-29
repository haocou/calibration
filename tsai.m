close all
clear
clc

[HA, HB] = LoadData();
%使用无噪声的数据求解手眼矩阵
[~,~,num]=size(HB);
Hce = HandEyeTsai(HA, HB,num);

%构造数据
function [HA, HB] = LoadData()

load('robot_tfs.mat');
load('pattern_tfs.mat');

%how to choose?????????????????
numI = [1,  88,  10,  66,  32, 40, 69, 100,  87,  93,  107, 145];
numJ = [88, 36, 119, 133, 93, 77, 198, 179, 123, 182, 147, 90];

[~, col] = size(numI);
HA=zeros(4,4,col);
HB=zeros(4,4,col);

for i=1:col
    HA(1:4,1:4,i) = reshape(A(numJ(i),1:4,1:4), 4,4) \ reshape(A(numI(i),1:4,1:4), 4,4);
    HB(1:4,1:4,i) = reshape(pattern_tfs(numJ(i),1:4,1:4), 4,4) / reshape(pattern_tfs(numI(i),1:4,1:4), 4,4);
end
end

function Hx = HandEyeTsai(HA, HB, num)
%num为使用输入数据中多少组数据进行标定

%Hx用于返回手眼标定结果
Hx = eye(4,4);

%根据Tsai方法生成最小二乘法求解的方程Ax = B，求解手眼矩阵的旋转矩阵
A = zeros(num*3, 3);
B = zeros(num*3, 1);

for i = 1:num
   %使用罗德里格斯计算旋转矩阵对应的旋转轴
%     rgij = rodrigues( HA(1:3, 1:3, i) );        
%     rcij = rodrigues( HB(1:3, 1:3, i) ); 
    rgij = rotationMatrixToVector( HA(1:3, 1:3, i) )';        
    rcij = rotationMatrixToVector( HB(1:3, 1:3, i) )'; 

    theta_gij = norm(rgij);
    theta_cij = norm(rcij);                        
    rngij = rgij/theta_gij;
    rncij = rcij/theta_cij;
    
    %Tsai方法里定义的旋转轴
    Pgij = 2*sin(theta_gij/2)*rngij;
    Pcij = 2*sin(theta_cij/2)*rncij; 

   %构造求解手眼旋转矩阵对应旋转轴的求解方程
    S = Skew( Pgij + Pcij );
    b = Pcij - Pgij;
    
    A(3*i-2 : 3*i, 1:3) = S;
    B(3*i-2 : 3*i, 1) = b;
end
 
    %最小二乘法求解旋转矩阵对应的旋转轴，然后根据Tsai的定义计算旋转矩阵
    Pce1 = inv(A'*A)*A'*B;
    Pce   =  2*Pce1/sqrt( 1 + norm(Pce1) * norm(Pce1) ) ;
    I = eye(3,3);
    np2 = norm(Pce) * norm(Pce);
    %得到手眼关系的旋转矩阵
    Rx = ( 1 - np2*0.5 )*I + 0.5*( Pce*Pce' + sqrt(4 - np2)* Skew(Pce) );
  
    %构造求解手眼关系平移量的方程Ax = B
    A = zeros(num*3, 3);
    B = zeros(num*3, 1);
for i = 1:num
    RA =  HA(1:3, 1:3, i);
    TA =  HA(1:3, 4, i);
   
    TB =  HB(1:3, 4,i);
   
    A(3*i-2 : 3*i, 1:3) = RA - I;
    B(3*i-2 : 3*i, 1) = Rx * TB - TA;
end
  %最小二乘法求解手眼关系的平移量
    Tx = inv(A'*A)*A'*B;
    Hx(1:3,1:3) = Rx;
    Hx(1:3, 4) = Tx;
end


function S = Skew( P )
    S = zeros(3,3);
    S(1,2) = -P(3);
    S(1,3) = P(2);
    
    S(2,1) = P(3);
    S(2,3) = -P(1);
    
     S(3,1) = -P(2);
    S(3,2) = P(1);
end

function M = EulerTomatrix(roll, pich, yaw)
%zyx
    m1 = [ 1.0,    0.0,            0.0;
               0.0,    cos(roll),    -sin(roll);
               0.0,    sin(roll),     cos(roll)];

    m2 = [ cos(pich),   0.0,   sin(pich);
               0.0,            1.0,      0.0;
              -sin(pich),     0.0,    cos(pich);];

    m3 = [ cos(yaw),   -sin(yaw),     0.0;
               sin(yaw),    cos(yaw),     0.0;
               0.0,                 0.0,         1.0];
    M = m3*m2*m1;
end

function	[out,dout]=rodrigues(in)

% RODRIGUES	Transform rotation matrix into rotation vector and viceversa.
%		
%		Sintax:  [OUT]=RODRIGUES(IN)
% 		If IN is a 3x3 rotation matrix then OUT is the
%		corresponding 3x1 rotation vector
% 		if IN is a rotation 3-vector then OUT is the 
%		corresponding 3x3 rotation matrix
%

%%
%%		Copyright (c) March 1993 -- Pietro Perona
%%		California Institute of Technology
%%

%% ALL CHECKED BY JEAN-YVES BOUGUET, October 1995.
%% FOR ALL JACOBIAN MATRICES !!! LOOK AT THE TEST AT THE END !!

%% BUG when norm(om)=pi fixed -- April 6th, 1997;
%% Jean-Yves Bouguet

%% Add projection of the 3x3 matrix onto the set of special ortogonal matrices SO(3) by SVD -- February 7th, 2003;
%% Jean-Yves Bouguet

% BUG FOR THE CASE norm(om)=pi fixed by Mike Burl on Feb 27, 2007


[m,n] = size(in);
%bigeps = 10e+4*eps;
bigeps = 10e+20*eps;

if ((m==1) & (n==3)) | ((m==3) & (n==1)) %% it is a rotation vector
    theta = norm(in);
    if theta < eps
        R = eye(3);

        %if nargout > 1,

        dRdin = [0 0 0;
            0 0 1;
            0 -1 0;
            0 0 -1;
            0 0 0;
            1 0 0;
            0 1 0;
            -1 0 0;
            0 0 0];

        %end;

    else
        if n==length(in)  in=in'; end; 	%% make it a column vec. if necess.

        %m3 = [in,theta]

        dm3din = [eye(3);in'/theta];

        omega = in/theta;

        %m2 = [omega;theta]

        dm2dm3 = [eye(3)/theta -in/theta^2; zeros(1,3) 1];

        alpha = cos(theta);
        beta = sin(theta);
        gamma = 1-cos(theta);
        omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        A = omega*omega';

        %m1 = [alpha;beta;gamma;omegav;A];

        dm1dm2 = zeros(21,4);
        dm1dm2(1,4) = -sin(theta);
        dm1dm2(2,4) = cos(theta);
        dm1dm2(3,4) = sin(theta);
        dm1dm2(4:12,1:3) = [0 0 0 0 0 1 0 -1 0;
            0 0 -1 0 0 0 1 0 0;
            0 1 0 -1 0 0 0 0 0]';

        w1 = omega(1);
        w2 = omega(2);
        w3 = omega(3);

        dm1dm2(13:21,1) = [2*w1;w2;w3;w2;0;0;w3;0;0];
        dm1dm2(13: 21,2) = [0;w1;0;w1;2*w2;w3;0;w3;0];
        dm1dm2(13:21,3) = [0;0;w1;0;0;w2;w1;w2;2*w3];

        R = eye(3)*alpha + omegav*beta + A*gamma;

        dRdm1 = zeros(9,21);

        dRdm1([1 5 9],1) = ones(3,1);
        dRdm1(:,2) = omegav(:);
        dRdm1(:,4:12) = beta*eye(9);
        dRdm1(:,3) = A(:);
        dRdm1(:,13:21) = gamma*eye(9);

        dRdin = dRdm1 * dm1dm2 * dm2dm3 * dm3din;


    end;
    out = R;
    dout = dRdin;

    %% it is prob. a rot matr.
elseif ((m==n) & (m==3) & (norm(in' * in - eye(3)) < bigeps)...
        & (abs(det(in)-1) < bigeps))
    R = in;

    % project the rotation matrix to SO(3);
    [U,S,V] = svd(R);
    R = U*V';

    tr = (trace(R)-1)/2;
    dtrdR = [1 0 0 0 1 0 0 0 1]/2;
    theta = real(acos(tr));


    if sin(theta) >= 1e-4,

        dthetadtr = -1/sqrt(1-tr^2);

        dthetadR = dthetadtr * dtrdR;
        % var1 = [vth;theta];
        vth = 1/(2*sin(theta));
        dvthdtheta = -vth*cos(theta)/sin(theta);
        dvar1dtheta = [dvthdtheta;1];

        dvar1dR =  dvar1dtheta * dthetadR;


        om1 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';

        dom1dR = [0 0 0 0 0 1 0 -1 0;
            0 0 -1 0 0 0 1 0 0;
            0 1 0 -1 0 0 0 0 0];

        % var = [om1;vth;theta];
        dvardR = [dom1dR;dvar1dR];

        % var2 = [om;theta];
        om = vth*om1;
        domdvar = [vth*eye(3) om1 zeros(3,1)];
        dthetadvar = [0 0 0 0 1];
        dvar2dvar = [domdvar;dthetadvar];


        out = om*theta;
        domegadvar2 = [theta*eye(3) om];

        dout = domegadvar2 * dvar2dvar * dvardR;


    else
        if tr > 0; 			% case norm(om)=0;

            out = [0 0 0]';

            dout = [0 0 0 0 0 1/2 0 -1/2 0;
                0 0 -1/2 0 0 0 1/2 0 0;
                0 1/2 0 -1/2 0 0 0 0 0];
        else

            % case norm(om)=pi;
            if(0)

                %% fixed April 6th by Bouguet -- not working in all cases!
                out = theta * (sqrt((diag(R)+1)/2).*[1;2*(R(1,2:3)>=0)'-1]);
                %keyboard;

            else

                % Solution by Mike Burl on Feb 27, 2007
                % This is a better way to determine the signs of the
                % entries of the rotation vector using a hash table on all
                % the combinations of signs of a pairs of products (in the
                % rotation matrix)

                % Define hashvec and Smat
                hashvec = [0; -1; -3; -9; 9; 3; 1; 13; 5; -7; -11];
                Smat = [1,1,1; 1,0,-1; 0,1,-1; 1,-1,0; 1,1,0; 0,1,1; 1,0,1; 1,1,1; 1,1,-1;
                    1,-1,-1; 1,-1,1];

                M = (R+eye(3,3))/2;
                uabs = sqrt(M(1,1));
                vabs = sqrt(M(2,2));
                wabs = sqrt(M(3,3));

                mvec = ([M(1,2), M(2,3), M(1,3)] + [M(2,1), M(3,2), M(3,1)])/2;
                syn  = ((mvec > eps) - (mvec < -eps)); % robust sign() function
                hash = syn * [9; 3; 1];
                idx = find(hash == hashvec);
                svec = Smat(idx,:)';

                out = theta * [uabs; vabs; wabs] .* svec;

            end

            if nargout > 1
                fprintf(1,'WARNING!!!! Jacobian domdR undefined!!!\n');
                dout = NaN*ones(3,9);
            end
        end
    end

else
    error('Neither a rotation matrix nor a rotation vector were provided');
end
end
