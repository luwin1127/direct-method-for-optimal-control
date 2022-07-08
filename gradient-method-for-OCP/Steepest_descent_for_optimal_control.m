%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 说明：应用最速下降法求解最优控制问题
% 例子：《最优化与最优控制》 pp. 257 例13.1 
% 类型：无控制约束的最优控制问题
% 作者：李凌炜
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;close all;

%% 主函数部分
eps = 1e-3;
options = odeset('RelTol', 1e-4, 'AbsTol',1e-4);
t0 = 0; tf = 1;
t_segment = 50;
Tu = linspace(t0, tf, t_segment);	% 将时间离散化
u = -0.5*ones(1, t_segment);        % 控制量的初值猜测
initx = 10;                         % 状态变量初值
initp = 0.1;                          % 协态变量初值
max_iteration = 50;                 % 最大迭代次数
H_Norm = zeros(max_iteration,1);
step = 0.55*ones(max_iteration,1);
J = zeros(max_iteration,1);
d = zeros(max_iteration,max_iteration);

%% 算法开始
tic;
% 记录计算时长
% mc = 10;
% time_record = zeros(mc,1);
% 
% for index = 1:mc
for i = 1:max_iteration   
    disp('---------------------------');
    fprintf('%d iterarion times \n',i)
    % 1) 正向积分状态方程得到 x(t) 
    [Tx, X] = ode45(@(t,x) stateEq(t,x,u,Tu), [t0 tf], initx, options);
    
    % 2) 反向积分协态方程得到 p(t)
    x = X;
    [Tp, P] = ode45(@(t,p) costateEq(t,p,x,Tx), [tf, t0], initp, options);
    p = P;
    
    % 因为求得的协态变量按照时间反序排列，为保证 x 和 p 在同一时间序列下
    % 需要调整其顺序
    p =interp1(Tp,p,Tx);
    
    % 计算 dH/du，即目标泛函的梯度deltaH
    dH = gradient(p,Tx,u,Tu);
    H_Norm(i,1) = dH'*dH;
    
    % 计算目标函数的值
    J(i,1) = fun(tf,x,Tx,u,Tu);
    fprintf('The value of the objective function is %.4f \n',J(i,1))
    
    % 如果梯度 dH/du < epsilon，则退出循环
    if H_Norm(i,1) < eps
        % 显示最后一次的目标函数值
%         fprintf('The value of the objective function is %.4f \n',J(i,1))
        break;
    else
        % 否则，更新控制量，进入下一次迭代
        u_old = u;
        % 计算更新步长
        % 我发现步长的影响更大，使不使用armijo线搜索来搜索步长对结果影响不大（beta = 0.55）
        % step = 0.1*ones(max_iteration,1) 时，用梯度法更新 u 的迭代次数为 39 次，共轭梯度法为 14 次
        % step = 0.55*ones(max_iteration,1) 时，用梯度法更新 u 的迭代次数为 5 次，共轭梯度法为 6 次
        % 而step = 0*ones(max_iteration,1) 时，用armijo线搜索来搜索步长，
        % 则梯度法更新 u 的迭代次数为 6 次，共轭梯度法为 13 次
        step(i,1) = armijo(u_old, -dH, tf, x, Tx, Tu);
        u = control_update_by_gradient(dH,Tx,u_old,Tu,step(i,1));
%         [u,d] = control_update_by_conjugate_gradient(H_Norm,dH,d,Tx,u_old,Tu,step(i,1),i);
    end
end

if i == max_iteration
    disp('---------------------------');
    disp('已达最大迭代次数，停止算法.');
end
toc;
% time_record(index) = toc;
% end
% time_cal = sum(time_record)/mc;
%% 画图
window_width = 500;
window_height = 416;
figure('color',[1 1 1],'position',[300,200,window_width,window_height]);
plot(Tx, X, 'x-','LineWidth',1.5);hold on;
plot(Tu, u, 'x-','LineWidth',1.5);
xlabel('Time');
ylabel('State & control');
set(gca,'FontSize',15,...
        'FontName','Times New Roman',...
        'LineWidth',1.5);
legend('$x(t)$','$u(t)$',...
        'Location','Northeast',...
        'FontSize',10,...
        'interpreter','latex');

count = find(J==0);
if isempty(count)
    t_plot = (1:length(J))';
    J_plot = J;
else
    if J(count(1)-1)-J(count(1)) > 1 
        t_plot = (1:count(1)-1)';
        J_plot = J(1:count-1);
    else
        t_plot = (1:length(J))';
        J_plot = J;
    end
end
figure('color',[1 1 1],'position',[300+window_width,200,window_width,window_height]);
plot(t_plot,J_plot,'x-','LineWidth',1.5);
xlabel('Iterations');
ylabel('J');
set(gca,'FontSize',15,...
        'FontName','Times New Roman',...
        'LineWidth',1.5);
legend('$J$',...
       'Location','Northeast',...
       'FontSize',10,...
       'interpreter','latex');

%% 子函数定义部分
% 状态方程
function dx = stateEq(t,x,u,Tu)
    u = interp1(Tu,u,t);
    dx = -x^2+u;
end

% 协态方程
function dp = costateEq(t,p,x,Tx)
    x = interp1(Tx,x,t);
    dp = -x + 2*p.*x;
end

% 用梯度法求解 dH/du 的表达式
function dH = gradient(p,Tx,u,Tu)
    u = interp1(Tu,u,Tx);
    dH = p + u;
end

% 更新控制量（梯度法）
function u_new = control_update_by_gradient(dH,tx,u,tu,step)
    beta = 0.55;
    dH = interp1(tx,dH,tu);
    % 用armijo准则搜索步长
%     u_new = u - beta^step*dH;
    % 不用armijo准则搜索步长
    u_new = u - step*dH;
end

% 更新控制量（共轭梯度法）
function [u_new,d] = control_update_by_conjugate_gradient(H_Norm,dH,d,tx,u,tu,step,i)
    beta = 0.55;
    dH = interp1(tx,dH,tu);
    if i == 1
        d(:,i) = -dH;
    else
        beta = H_Norm(i)/H_Norm(i-1);
        d(:,i) = -dH' + beta*d(:,i-1);
    end
    % 用armijo准则搜索步长
    u_new = u + beta^step*d(:,i)';
    % 不用armijo准则搜索步长
%     u_new = u + step*d(:,i)';
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 说明：armijo准则非精确线搜索程序，用来求解迭代步长
% 文献：《最优化方法及其Matlab程序设计》 pp. 26
% 作者：袁亚湘
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mk = armijo(xk, dk, tf, x, Tx, Tu)
dk = -interp1(Tx,dk,Tu);
beta = 0.55;
sigma = 0.5;
m = 0;
mmax = 20;
while m <= mmax
    if fun(tf,x,Tx,xk+beta^m*dk,Tu) <= fun(tf,x,Tx,xk+beta^m*dk,Tu)+sigma*beta^m*gfun(tf,xk,Tu)*dk'
        m = m+1;
        break;
    end
    m = m+1;
end
mk = m;
% 这个 end 归属于 function 
end

function J = fun(tf,x,Tx,u,Tu)
    % 2022-07-01 的代码
    % 写打靶法的时候发现梯度法求目标函数的方法错了
    % 更新了求目标函数的方式，我认为这个求法才是对的
    u = interp1(Tu,u,Tx);
    L = x.^2 / 2 + u.^2 /2;
    J = tf*trapz(Tx,L);
    % 2022-06-25 的代码
%     J = tf*(((x')*x)/length(Tx) + (u*u')/length(Tu));
end

function gJ = gfun(tf,u,Tu)
    gJ = tf*(2*u)/length(Tu);
end