%--------------------------------------------------------------------------
% This code demonstrates an example of solving constrained optimization problem 
% with multiple shooting method.
% Author: Vinh Quang Nguyen - University of Massachusetts, Amherst
%--------------------------------------------------------------------------
% 说明：应用多重打靶法求解无约束最优控制问题
% 例子：《最优化与最优控制》 pp. 257 例13.1 
% 类型：无控制约束的最优控制问题
% 时间：2022/07/05
%--------------------------------------------------------------------------
clear;clc;close all;

%% 01 初始参数设置
p.ns = 1; p.nu = 1;                     % 状态量个数和控制量个数
p.t0 = 0; p.tf = 1;                     % 初始时间和终止时间
p.x0 = 10;                              % 初始条件

% 多重打靶法参数设置
p.N = 20;                               % 打靶点数 => (N-1) 个子时间区段
p.M = 4;                                % 每个子时间区段包含的打靶点
p.t = linspace(p.t0,p.tf,p.N);          % 时间序列

% 设置状态量和控制量的索引
p.x_index = 1:p.N;
p.u_index = p.N+1:2*p.N;
%% 02 求解算法
% 设置初值
y0 = ones((p.ns + p.nu)* p.N, 1);
% 设定求解器设置
options = optimoptions('fmincon','Display','Iter','Algorithm','sqp','MaxFunEvals',1e5); 
% mc = 10;
% time_record = zeros(mc,1);
% 
% for index = 1:mc
tic;
[X,fval,exitflag,output] = fmincon(@(y) objfun(y, p),y0,[],[],[],[],[],[],@(y) noncon(y, p),options);
toc;
% time_record(index) = toc; 
% end
% time_cal = sum(time_record)/mc;

%% 03 处理数据
p.x = reshape(X(p.x_index), [], p.ns);
p.u = reshape(X(p.u_index), [], p.nu);
% 尝试外推最后一个时间点的值
% p.u(p.N) = (p.t(p.N)-p.t(p.N-1)) * (p.u(p.N-1)-p.u(p.N-2)) / (p.t(p.N-1)-p.t(p.N-2)) + p.u(p.N-1);
% 尝试平滑第一个时间点的值
u0 = p.u(1+1) + (p.u(1+1)-p.u(2+1))*(p.t(0+1)-p.t(1+1))/(p.t(1+1)-p.t(2+1));
p.u(1) = u0;
%% 04 画图
window_width = 500;
window_height = 416;

% 状态量和控制量
k = 1;
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
% save(['.\','multi_shooting_method.mat']);
%% 子函数  
% 目标函数
function f = objfun(y,p)
    % 得到状态量和控制量
    x = y(p.x_index);
    u = y(p.u_index);
    % 为了保证 x 和 u 的维度一致，需要对 u 进行外推，得到 u 在末端时刻的值
    % 利用简单的线性外推方法进行外推
    % k = (x2-x1)/(t2-t1) == (x3-x2)/(t3-t2)
    % x3 = (t3-t2)(x2-x1)/(t2-t1) + x2
%     N = p.N;
%     t = p.t;
%     u(N) = (t(N)-t(N-1)) * (u(N-1)-u(N-2)) / (t(N-1)-t(N-2)) + u(N-1);
    L = u.^2/2 + x.^2/2;            % 积分项
    f = trapz(p.t,L);               % 计算目标函数
end

% 状态方程
function dy = state_eq(y,u)
    dy = -y^2 + u;
end

% 约束条件
function [c,ceq] = noncon(y,p)
    % 得到状态量和控制量
    x = reshape(y(p.x_index),[],p.ns);
    u = reshape(y(p.u_index),[],p.nu);
    
    % 时间步长
    h = p.tf/(p.N-1)/(p.M-1);
    
    % 每次子时间区段进行单次打靶法
    states_at_nodes = zeros(p.N, p.ns);
    for i = 1:p.N-1
       x0 = x(i,:);
       u0 = u(i,:);
       states = zeros(p.M,p.ns);
       states(1,:) = x0;
       for j =1:p.M-1
           k1 = state_eq(states(j,:), u0);
           k2 = state_eq(states(j,:) + h/2 * k1, u0);
           k3 = state_eq(states(j,:) + h/2 * k2, u0);
           k4 = state_eq(states(j,:) + h * k3, u0);
           states(j+1,:) = states(j,:) + h/6*(k1 + 2*k2 + 2*k3 + k4);
       end
       states_at_nodes(i+1,:) = states(end,:);
    end
    
    % 保证各区段起始点的连续性
    ceq_temp = x(2:end,:) - states_at_nodes(2:end,:);
    
    % 把初始时刻的状态约束放到 ceq 中
    ceq_temp = [ceq_temp; x(1,:) - p.x0];
    ceq = reshape(ceq_temp, [], 1);
    
    % 不等式约束
    c = [];
end