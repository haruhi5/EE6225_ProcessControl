%% EE6225 process control assignment 1
% Rin Misaka
%% system response
sys = tf(1.5,[1 11 45 85 78 36],'outputdelay',2) % exp(-2s) - time delay
t = -1:0.1:20;
response = step(sys,t);
unitstep = t>=0;
figure
plot(t,unitstep); grid on;
title('Unit Step Response');

% 数组补0
response = [zeros(length(t)-length(response),1); response];
len = length(response);

figure
plot(t,response);
xlabel('time'); ylabel('response');
title('Step Response');
hold on; grid on;
%plot(t,unitstep,'r--');
%steady-value
for i=32:len-1
    if(abs(response(i)-response(i-1))<0.0000001 && abs(response(i)-response(i+1))<0.0000001)
        steady_index = i;
        break;
    end
end
%steady_value = sum(response(len-40:len))/length(response(len-40:len)); %y_∞
steady_value = sum(response(steady_index:len))/length(response(steady_index:len)); %y_∞

plot(t,steady_value*ones(1,length(t)),'r--');

%% two points method
fprintf("Two points method\n")
K_gain = steady_value/1;
point1 = 0.284*steady_value;
point2 = 0.632*steady_value;
fprintf("Point 1 (28.4%%): %f\nPoint 2 (63.2%%): %f\n",point1,point2);
%找到对应的两个点的坐标
label1 = 0; label2 = 0; temp = 0;
for i=1:len
    if(response(i)-point1<0.001)
        if(temp<i) 
            label1 = temp;
        end
        temp=i;
    end
    if(response(i)-point2<0.001)
        label2 = i;
    end
end
t1 = t(label1);t2 = t(label2);
% 求解T，L
T_towpoints = 1.5*(t2-t1);
L_twopoints = 0.5*(3*t1-t2);
fprintf("Two points method‘s solution:\nT(time constant): %f\nL(delay): %f\n\n", T_towpoints,L_twopoints);

% plot the predict two point method model
sys_twopoint = tf(K_gain,[T_towpoints 1],'outputdelay',L_twopoints)
response_twopoint = step(sys_twopoint,t);
response_twopoint = [zeros(length(t)-length(response_twopoint),1); response_twopoint];

figure
plot(t,response);hold on; grid on;
plot(t,response_twopoint);
plot(t,steady_value*ones(1,length(t)),'r--');
xlabel('time'); ylabel('response');
title('Step Response of two point method');
figure
nyquist(sys,sys_twopoint,{0.0001,pi/4});


%% log method
fprintf("Log method\n")
K_gain = steady_value/1;
%ln[(y_∞-y)/y_∞]=L/T - t/T
%y_temp = ln[(y_∞-y)/y_∞]
temp = (steady_value-response)/steady_value;
%choose some points to fit line % 有些是负数？delay的时间里响应为0
fit_point = temp(label1:label2);
y_temp = log(fit_point);
figure
plot(t(label1:label2),y_temp,'o'); 
grid on; hold on;
% for i=1:len
%    if(temp(i)<0) response_positiveindex = i-1; break;
%    end
% end y_temp = temp(1:response_positiveindex); y_temp = log( y_temp ); figure plot(t(1:length(y_temp)),y_temp,'o')

%fit
t_fit = t(label1:label2)';
%b_fit = [ones(length(y_temp),1) y_temp]\t_fit; 
% 除反了 线性拟合: y = β0+β1*t -> 变为矩阵形式 !!!
b_fit = [ones(length(t_fit),1) t_fit]\y_temp;
y_fit = b_fit(2)*(t)+b_fit(1);
plot(t,y_fit);
L_logmethod = -b_fit(1)/b_fit(2);
T1_logmethod = L_logmethod/b_fit(1);%T2_logmethod = -1/b_fit(2);
fprintf("Log method‘s solution:\nA(gain):%f\nT(time constant): %f\nL(delay): %f\n\n",K_gain,T1_logmethod, L_logmethod);

% plot predict log method model
sys_logmethod = tf(K_gain,[T1_logmethod 1],'outputdelay',L_logmethod)
response_logmethod = step(sys_logmethod,t);
response_logmethod = [zeros(length(t)-length(response_logmethod),1); response_logmethod];

figure
plot(t,response);hold on; grid on;
plot(t,response_logmethod);
plot(t,steady_value*ones(1,length(t)),'r--');
xlabel('time'); ylabel('response');
title('Step Response of log method');
figure
nyquist(sys,sys_logmethod,{0.0001,pi/4});



%% area method
fprintf("Area method\n")
y_infinit_area = t(steady_index)*steady_value;
%yt_area = trapz(0.1,response(1:71))+(14.9-6)*steady_value; % A1 area 这不是A1！！！
yt_area = trapz(0.1,response(1:160)); % A1 area 这不是A1！！！A1是积分到Tar
Tar = (y_infinit_area-yt_area)/K_gain; % Tar = T+L
A1 = trapz(0.1,response(1:53))
T_areamethod = exp(1)*A1/K_gain;
L_areamethod = Tar-T_areamethod;
fprintf("Area method‘s solution:\nA(gain):%f\nT(time constant): %f\nL(delay): %f\n\n",K_gain,T_areamethod, L_areamethod);

% plot the predict area method model
% plot predict log method model
sys_areamethod = tf(K_gain,[T_areamethod 1],'outputdelay',L_areamethod)
response_areamethod = step(sys_areamethod,t);
response_areamethod = [zeros(length(t)-length(response_areamethod),1); response_areamethod];

figure
plot(t,response);hold on; grid on;
plot(t,response_areamethod);
plot(t,steady_value*ones(1,length(t)),'r--');
xlabel('time'); ylabel('response');
title('Step Response of area method');
figure
nyquist(sys,sys_areamethod,{0.0001,pi/4});
