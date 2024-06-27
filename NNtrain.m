clc,clear
%% Import data
load('dataset.mat');
test_data_num = randperm(size(dataset,1)); 
test_data_num = test_data_num(1:round(size(dataset,1)/10))';
train_data_num = setdiff([1:size(dataset,1)]',test_data_num);
object_num = 6; inp = dataset(:,1:5); outp = dataset(:,object_num);
%% Data normalisation
[inp_i,ps_i] = mapminmax(inp(train_data_num,:)');
[outp_t,ps_t] = mapminmax(outp(train_data_num,:)');
x = inp_i; t = outp_t;
% save('ps_i','ps_i'); save('ps_t','ps_t'); 
%% Choose a Training Function
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.
hiddenLayerSize = [40 40 10]; % Create a Fitting Network
net = fitnet(hiddenLayerSize,trainFcn);
net.layers{1}.transferFcn = 'tansig';
net.layers{2}.transferFcn = 'tansig';
net.layers{3}.transferFcn = 'tansig';
%% Setup Division of Data for Training, Validation, Testing elliotsig
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 20/100;
net.divideParam.testRatio = 10/100;
%% Train the Network
[net,tr] = train(net,x,t);
% save('net','net'); 
%% Test the Network
x_test = mapminmax('apply',inp(test_data_num,:)',ps_i);
y_test = net(x_test);
yout = mapminmax('reverse',y_test,ps_t);
e = gsubtract(outp(test_data_num,:)',yout);
rmse = sqrt((1/size(e,2))*sum(e.^2));
% View the Network
% view(net)
% Uncomment these lines to enable various plots.
figure(1)
X = outp(test_data_num,:); Y = yout;
plot([X],[Y],'o','color',[121 42 5]./255);
hold on
plot([Y],[Y],'linewidth',2,'color',[102 147 72]./255);
hold off
set(gca,'YLim',[min(Y) max(Y)]); set(gca,'XLim',[min(Y) max(Y)]);
set(gca,'YTick',[min(Y) max(Y)]); set(gca,'XTick',[min(Y) max(Y)]); set(gca,'FontSize',20);
h1 = legend('Data','Fit');
set(h1,'Location','NorthWest','Box','off','FontSize',15);
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%  Please feel free to contact us with any questions! 
%  - Chuang Ma, Dalian University of Technology
%  - Yichao Zhu, Dalian University of Technology
%  - Xu Guo, Dalian University of Technology
%  - chuangma@mail.dlut.edu.cn / yichaozhu@dlut.edu.cn / guoxu@dlut.edu.cn
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%