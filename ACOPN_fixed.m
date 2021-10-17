close all
clc
clear
%%初始化标识和变迁节点%%
nodes_data=cell(0);
M_1=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; %网中25个标识对应的符号
M_2=[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_3=[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_4=[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_5=[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_6=[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_7=[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_8=[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_9=[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_10=[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_11=[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_12=[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0];
M_13=[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0];
M_14=[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
M_15=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
M_16=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0];
M_17=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0];
M_18=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0];
M_19=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0];
M_20=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0];
M_21=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0];
M_22=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0];
M_23=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0];
M_24=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0];
M_25=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];
nodes_data(1,:)={1,[2,3],[1,1],[M_2;M_3]};  %根据变迁的邻近变迁的发生权构造变迁元胞组
nodes_data(2,:)={2,[4,5],[1,sqrt(2)],[M_4;M_5]};
nodes_data(3,:)={3,[6],[1],[M_6]};
nodes_data(4,:)={4,[7,8],[sqrt(2),1],[M_7;M_5]};
nodes_data(5,:)={5,[9,10],[sqrt(2),1],[M_8;M_7]};
nodes_data(6,:)={6,[11,12],[1,sqrt(2)],[M_9;M_10]};
nodes_data(7,:)={7,[13,14],[1,1],[M_8;M_13]};
nodes_data(8,:)={8,[9,10],[sqrt(2),1],[M_8;M_7]};
nodes_data(9,:)={9,[15],[1],[M_18]};
nodes_data(10,:)={10,[13,14],[1,1],[M_8;M_13]};
nodes_data(11,:)={11,[16,17],[1,1],[M_10;M_11]};
nodes_data(12,:)={12,[18],[1],[M_12]};
nodes_data(13,:)={13,[15],[1],[M_18]};
nodes_data(14,:)={14,[23,24],[sqrt(2),1],[M_14;M_15]};
nodes_data(15,:)={15,[29],[1],[M_19]};
nodes_data(16,:)={16,[18],[1],[M_12]};
nodes_data(17,:)={17,[19,20],[1,sqrt(2)],[M_16;M_17]};
nodes_data(18,:)={18,[21,22],[1,sqrt(2)],[M_18;M_19]};
nodes_data(19,:)={19,[27],[1],[M_17]};
nodes_data(20,:)={20,[28],[1],[M_20]};
nodes_data(21,:)={21,[29],[1],[M_19]};
nodes_data(22,:)={22,[30,31],[1,sqrt(2)],[M_23;M_24]};
nodes_data(23,:)={23,[25],[1],[M_21]};
nodes_data(24,:)={24,[26],[1],[M_14]};
nodes_data(25,:)={25,[33],[1],[M_22]};
nodes_data(26,:)={26,[25],[1],[M_21]};
nodes_data(27,:)={27,[28],[1],[M_20]};
nodes_data(28,:)={28,[32],[1],[M_23]};
nodes_data(29,:)={29,[30,31],[1,sqrt(2)],[M_23;M_24]};
nodes_data(30,:)={30,[35],[1],[M_24]};
nodes_data(31,:)={31,[36],[1],[M_25]};
nodes_data(32,:)={32,[35],[1],[M_24]};
nodes_data(33,:)={33,[34],[1],[M_25]};
nodes_data(35,:)={35,[36],[1],[M_25]};


%% 始末节点%%
node_start=1;
node_end=[34,36];
%%% 蚁群定义%%%%%
m=50;                               % 蚂蚁数量
n=size(nodes_data,1);               % 节点数量
alpha=1;                            % 信息素重要程度因子
beta=5;                             % 启发函数重要程度因子
Rho=0.5;                            % 信息素挥发因子
Q=1;                              % 信息素增加强度系数
%%迭代过程初始化定义%%%%
iter=1;                             % 迭代次数初值
iter_max=100;                       % 最大迭代次数
Route_best=cell(iter_max,1);        % 各代最佳路径
Length_best=zeros(iter_max,1);      % 各代最佳路径长度
Length_ave=zeros(iter_max,1);       % 各代路径平均长度
Place_best=cell(iter_max,1);        % 各代最佳路径访问的库所
%%将信息素、挥发因子一并放入nodes_data中%%%%%
Delta_Tau_initial=nodes_data(:,1:2);
for i=1:size(nodes_data,1)
    nodes_data{i,5}=ones(1,length(nodes_data{i,3}));         % 初始信息素均设置为1
    nodes_data{i,6}=1./nodes_data{i,3};                      % 启发函数设置为距离的倒数
    Delta_Tau_initial{i,3}=zeros(1,length(nodes_data{i,3})); % 信息素变化量均为0
end
%% 迭代寻找最佳路径%%%
while iter<iter_max
    route=cell(0);
    place=cell(0);
    for i=1:m                                               % 逐个蚂蚁进路径选择
        neighbor_allow=cell(0);
        node_step=node_start;
        path=node_step;
        path_M = 0;
        if node_step==node_start
            marking=M_1;
        else
            marking=nodes_data{node_step,4};
        end
        dist=0;
        while ~isequal(marking(end,:), M_25)
            neighbor=nodes_data{node_step,2};                    %寻找邻近节点
            neighbor_allow = [];
            for k=1:length(neighbor)
                if ~ismember(neighbor(k),path)
                    neighbor_allow(end+1) = neighbor(k);
                end
            end           
            if isempty(neighbor_allow)
                neighbor_allow=cell(0);
                node_step=node_start;
                path=node_step;
                if node_step==node_start
                    marking=M_1;
                else
                    marking=nodes_data{node_step,4};
                end
                dist=0;
                continue
            end
            P=neighbor_allow;                                  %计算下一个节点的访问概率
            for k=1:length(neighbor_allow)
                idx = find(neighbor_allow(k) == nodes_data{node_step,2});
                P(2,k)=nodes_data{node_step,5}(idx)^alpha*nodes_data{node_step,6}(idx)^beta;
            end
            P(2,:)=P(2,:)/sum(P(2,:));
            %%轮盘赌法选择下一个访问节点%%%%
            Pc=cumsum(P(2,:));
            Pc=[0,Pc];
            randnum=rand;
            for k=1:length(Pc)-1
                if randnum>Pc(k)&&randnum<Pc(k+1)
                    target_node=neighbor_allow(k);
                end
            end
            %%%计算单步距离%%%
            idx=find(nodes_data{node_step,2}==target_node);
            dist=dist+nodes_data{node_step,3}(idx);
            %%%更新路径、节点以及标识%%%
            path(end+1)=target_node;                         % 更新路径集合
            marking(end+1,:)=nodes_data{node_step,4}(idx,:); % 更新标识
            node_step=target_node;                           % 更新下一个目标节点变迁
        end
        Length(i,1)=dist;                                    % 存放第i只蚂蚁的累计距离和对应路径
        route{i,1}=path;
        place{i,1}=marking;
    end
    %%计算这一代的m只蚂蚁中最短距离和对应路径%%%
    if iter==1
        [min_Length,min_index]=min(Length);
        Length_best(iter)=min_Length;
        Length_ave(iter)=mean(Length);
        Route_best{iter,1}=route{min_index,1};
        Place_best(iter,1)=place(min_index,1);
    else
        [min_Length,min_index]=min(Length);
        Length_best(iter)=min(Length_best(iter-1),min_Length);
        Length_ave(iter)=mean(Length);
        if Length_best(iter)==min_Length
            Route_best{iter,1}=route{min_index,1};
            Place_best(iter,1)=place(min_index,1);
        else
            Route_best{iter,1}=Route_best{iter-1,1};
            Place_best(iter,1)=Place_best(iter-1,1);
        end
    end
    %%%更新信息素%%%
    Delta_Tau=Delta_Tau_initial;
    for i=1:m                                                   % 逐个蚂蚁计算
        for j=1:length(route{i,1})-1
            node_start_temp=route{i,1}(j);
            node_end_temp=route{i,1}(j+1);
            idx=find(Delta_Tau{node_start_temp,2}==node_end_temp);
            if ismember(route{i,1}(j),[1,2,3,5,6,9,12,15,18,22,29,31,35]);
                Delta_Tau{node_start_temp,3}(idx)=Delta_Tau{node_start_temp,3}(idx)+1/8.2426;
            else
                Delta_Tau{node_start_temp,3}(idx)=Delta_Tau{node_start_temp,3}(idx);
            end
        end
    end
    %%考虑挥发因子，更新信息素%%%
    for i=1:size(nodes_data,1)
        nodes_data{i,5}=(1-Rho)*nodes_data{i,5}+Rho*Delta_Tau{i,3};
    end
    iter=iter+1;
end
% %绘图、结果%%%
figure
plot(1:iter_max,Length_best,'b',1:iter_max,Length_ave,'r');
legend('最短距离','平均距离');
xlabel('迭代次数');
ylabel('距离');
title('各代最短距离与平均距离对比');
%% 最优路径%%%
[dist_min,idx]=min(Length_best(1:end-1));
path_opt=Route_best{idx,1};
marking_opt=Place_best{idx,1};
%%将marking_opt转为字符输出%%%
M_set = cell(0);
M_seq = [];
for i = 1:size(marking_opt,1)
    idx = find(marking_opt(i,:) == 1);
    M_set{i,1} = strcat('M_',num2str(idx));
    M_seq = [M_seq, strcat(M_set{i,1},',')];
end
disp('机器人需要走的最短路径长度为：')
disp(dist_min)
disp('机器人经过的最短路径库所顺序为：')
disp(M_seq(1:end-1))
disp('机器人经过的最短路径变迁序列为：')
disp(path_opt)