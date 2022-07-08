%--------------------------------------------------------------------------
% Method_Pseudospectral.m
% Attempt to solve the Bryson-Denham problem using a pseudospectral method
% (namely LGL nodes and Gaussian quadrature) 
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 说明：应用伪谱法求解无约束最优控制问题
% 例子：《最优化与最优控制》 pp. 257 例13.1 
% 类型：无控制约束的最优控制问题
% 时间：2022/07/07
%--------------------------------------------------------------------------
clear;clc;close all;

%% 01 初始参数设置
p.ns = 1; p.nu = 1;                     % 状态量个数和控制量个数
p.t0 = 0; p.tf = 1;                     % 初始时间和终止时间
p.x0 = 10;                              % 初始条件

p.nt = 50;                              % 配置点个数
% p.method = 'LGL';                       % 选择配点方式：Legendre-Gauss-Lobatto(LGL)
p.method = 'CGL';                       % 或者 Chebyshev-Gauss-Lobatto(CGL)
% p.method = 'LGR';                       % 或者 Legendre-Gauss-Radau(LGR)
if strcmp(p.method,'LGL')
    p.tau = LGL_nodes(p.nt-1);          % 把时间区间离散到[-1,1]中
    p.D   = LGL_Dmatrix(p.tau);         % 微分矩阵
    p.w   = LGL_weights(p.tau);         % 配置点权重
elseif strcmp(p.method,'CGL')    
    p.tau = CGL_nodes(p.nt-1);          % 把时间区间离散到[-1,1]中
    p.D   = CGL_Dmatrix(p.tau);         % 微分矩阵
    p.w   = CGL_weights(p.tau);         % 配置点权重
elseif strcmp(p.method,'LGR')
    [p.tau,p.w] = LGR_nodes_and_weights(p.nt-1);
    p.D   = LGR_Dmatrix(p.tau);         % 微分矩阵
end

% 将变量离散为 x = [x,u]
p.x_index = 1:p.nt;
p.u_index = p.nt+1:2*p.nt;

%% 02 求解算法
% 给[x;u]赋初值
y0 = zeros(p.nt*(p.ns+p.nu),1);
% 求解器设置
options = optimoptions(@fmincon,'display','iter','MaxFunEvals',1e5,'algorithm','sqp');

mc = 10;
time_record = zeros(mc,1);

for index = 1:mc
tic;
[y,fval] = fmincon(@(y) objective(y,p),y0,[],[],[],[],[],[],@(y) constraints(y,p),options);
% toc;
time_record(index) = toc; 
end
time_cal = sum(time_record)/mc;

% 得到求解结果
p.x = y(p.x_index);                     % 状态量
p.u = y(p.u_index);                     % 控制量
p.t = tau2t(p);                         % 把时间变量从 tau 转化为 t 

%% 03 画图    
window_width = 500;
window_height = 416;

% 状态量
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
% save(['..\','pseudospectral_method.mat']);
%% 04 子函数部分
% 目标函数
function J = objective(y,p)
    % 得到状态量和控制量
    x = y(p.x_index);
    u = y(p.u_index);
    
    L = u.^2 + x.^2;                    % 积分项
    J = (p.tf-p.t0)/2 * dot(p.w,L)/2;   % 计算目标函数
end

% 约束函数
function [c,ceq] = constraints(y,p)
    % 得到状态量和控制量
    x = y(p.x_index);
    u = y(p.u_index);
    % 状态方程约束
    Y = x; F = (p.tf-p.t0)/2*(-x.^2+u);
    % 初始状态约束
    ceq1 = x(1) - p.x0;
    ceq2 = p.D*Y - F;
    % 输出约束条件    
    ceq = [ceq1;ceq2(:)];
    c = [];
end

% 把 tau 转换为 t
function t = tau2t(p)
    % 按照《天基对地打击武器轨道规划与制导技术》2013 pp.79 式(4.25)编写如下代码
    t = (p.tf-p.t0)*p.tau/2+(p.tf-p.t0)/2;
end