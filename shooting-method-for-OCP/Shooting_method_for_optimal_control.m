%--------------------------------------------------------------------------
% Method_SingleShooting.m
% Attempt to solve the Bryson-Denham problem using a single shooting method
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 说明：应用直接打靶法求解无约束最优控制问题
% 例子：《最优化与最优控制》 pp. 257 例13.1 
% 类型：无控制约束的最优控制问题
% 时间：2022/07/01
%--------------------------------------------------------------------------
clear;clc;close all;

%% 01 初始参数设置
p.ns = 1; p.nu = 1;                     % 状态量个数和控制量个数
p.t0 = 0; p.tf = 1;                     % 初始时间和终止时间
p.x0 = 10;                              % 初始条件

% 直接打靶法参数设置
p.nt = 20;                              % 打靶点参数设置
p.t = linspace(p.t0,p.tf,p.nt)';        % 时间区间

% 将控制量离散
p.u_index = 1:p.nt;

%% 02 求解算法
% 给控制量赋初值
p.u = -0.5*ones(p.nt,1);
u0 = p.u;                               % 控制量的初值猜测
% mc = 10;                                % 重复仿真次数
% time_record = zeros(mc,1);

% for index = 1:mc
% 求解控制量
options = optimoptions(@fmincon,'display','iter','MaxFunEvals',1e5,'algorithm','sqp');
tic;
[u,fval,exitflag,output] = fmincon(@(u) objective(u,p),u0,[],[],[],[],[],[],[],options);
toc;
p.u = u;
time_record = toc;
% time_record(index) = toc;
% end
% time_cal = sum(time_record)/mc;

% 再进行一次仿真得到数据
[~,Y] = ode45(@(t,y) derivative(t,y,p),p.t,p.x0);
p.x = Y(:,1);

%% 画图
k = 0;
window_width = 500;
window_height = 416;

% 状态量
figure('color',[1 1 1],'position',[300+k*window_width,300,window_width,window_height]);
plot(p.t, p.x, 'o-', 'LineWidth',1.5);hold on;
plot(p.t, p.u, 'x-', 'LineWidth',1.5);
xlabel('Time');
ylabel('State & control');
set(gca,'FontSize',15,...
        'FontName','Times New Roman',...
        'LineWidth',1.5);
legend('$x(t)$','$u(t)$',...
        'Location','Northeast',...
        'FontSize',10,...
        'interpreter','latex');

% 保存数据
% .\ 下一级文件夹
% ..\ 上一级文件夹
% save(['..\','shooting_method.mat']);

%% 子函数部分
% 目标函数
function f = objective(u_obj,p)
    p.u = u_obj(p.u_index);
    [~,Y] = ode45(@(t,y) derivative(t,y,p),p.t,p.x0); % 仿真得到时序状态量
    x = Y;                          % 状态量
    u = u_obj(p.u_index);           % 控制量
    L = u.^2/2 + x.^2/2;            % 积分项
    f = trapz(p.t,L);               % 计算目标函数
end

% 状态方程
function dy = derivative(t,y,p)
    % 使用 interp1qr() 进行插值
    u = interp1qr(p.t,p.u,t);
    % 使用 interp1q() 函数进行插值
    % 发现 interp1q() 的速度比 interp1() 快，和nterp1qr() 速度一样
%     u = interp1q(p.t,p.u,t);
    % 使用 interp1() 函数进行插值
%     u = interp1(p.t,p.u,t);    
    dy = -y^2 + u;
end
