%% EE6225 process control assignment 1-2-3
% Wu Tanghong
%% system response
sys = tf(1.5,[1 11 45 85 78 36],'outputdelay',2) % exp(-2s) - time delay
t = 0:0.1:20;
response = step(sys,t);
unitstep = t>=0;
% figure
% plot(t,unitstep); grid on;
% title('Unit Step Response');
figure
plot(t,response);
xlabel('time'); ylabel('response');
title('Step Response');
hold on; grid on;

len = length(response);

for i=32:len-1
    if(abs(response(i)-response(i-1))<0.0000001 && abs(response(i)-response(i+1))<0.0000001)
        steady_index = i;
        break;
    end
end
%steady_value = sum(response(len-40:len))/length(response(len-40:len)); %y_∞
steady_value = sum(response(steady_index:len))/length(response(steady_index:len)); %y_∞

plot(t,steady_value*ones(1,length(t)),'r--');

%% Open loop using LSM - time domain
fprintf("Open loop using LSM - time domain\n")

A = 1;
for i=1:len
   if(response(i)>0.0001)
       nonzero = i;
       break;
   end
end
phai = zeros(len-nonzero+1,3);
for i=1:len-nonzero
   % 要从t>L时开始！！！
   phai(i,:) = [-trapz(0.1,response(nonzero-1:nonzero-1+i)) -A A*t(nonzero-1+i)]; % 去掉0不影响
   %test(i,:) = [-trapz(0.1,response(22:i)) -A A*t(i)]; 
end

theta = (inv(phai'*phai))*phai'*response(nonzero:len);
a1 = theta(1);
b1 = theta(3);
L = theta(2)/theta(3);

sys_lsmtime = tf(b1,[1 a1],'outputdelay',L) % exp(-2s) - time delay 这个b1和a1位置对的
response_lsmtime = step(sys_lsmtime,t);
plot(t,response_lsmtime);
figure
nyquist(sys,sys_lsmtime,{0.0001,pi/4});
fprintf("Open loop using LSM - time domain method‘s solution:\nT(time constant): %f\nL(delay): %f\n\n", 1/a1,L);

% phai2 = zeros(len-22,3);
% 
% for i=1:len-22
%    %phai2(i,:) = [-(response(i)+response(i+1))*0.1/2 -A A*t(i)]; % 去掉0不影响
%    phai2(i,:) = [-trapz(0.1,response(21:i+21)) -A A*t(i+21)]; 
% end
% theta = (inv(phai2'*phai2))*phai2'*response(22:len);
% a1 = theta(1);
% b1 = theta(3);
% L = theta(2)/theta(3);
% 
% sys_lsmtime = tf(b1,[1 a1],'outputdelay',L) % exp(-2s) - time delay
% response_lsmtime = step(sys_lsmtime,t);
% plot(t,response_lsmtime);
%% Open loop using LSM - frequence domain
fprintf("Open loop using LSM - frequence domain\n")
A=1;
len = length(response);

I_sin = 0;
I_cos = 0;

syms w;
for j=1:len
        I_sin = I_sin+( response(j)-steady_value )*sin( w *t(j))*0.1; % 每个wi对应的积分
        I_cos = I_cos+( response(j)-steady_value )*cos( w *t(j))*0.1;
end

Gjw(w) = (steady_value+w*I_sin+1j*w*I_cos)/A;

% recursive
M=10; % PPT 3
Gjw_mag = zeros(M+1,1);
phai_g = zeros(M+1,1);
omiga = zeros(M+1,1);
omiga(1)=0;phai_g(1)=0;Gjw_mag(1)=abs(Gjw(omiga(1)));
omiga(2)=0.001;phai_g(2)=angle(Gjw(omiga(2)));Gjw_mag(2)=abs(Gjw(omiga(2)));
for i=3:M+1
omiga(i)=omiga(i-1)-((i-1)*pi/(M-1)+phai_g(i-1))*(omiga(i-1)-omiga(i-2))/(phai_g(i-1)-phai_g(i-2));

phai_g(i)=angle(Gjw(omiga(i)));
Gjw_mag(i)=abs(Gjw(omiga(i)));
end
%polar(phai_g,Gjw_mag);
phai=[-omiga.^2.*Gjw_mag.^2 ones(M+1,1)];

theta = (inv(phai'*phai))*phai'*(Gjw_mag.^2)
a1 = sqrt(theta(1))
b0 = sqrt(theta(2))

L = (inv(omiga'*omiga))*omiga'*(-phai_g-angle(1+1j*a1*omiga))

sys_lsmfreq = tf(b0,[a1 1],'outputdelay',L) 
response_lsmfreq = step(sys_lsmfreq,t);
plot(t,response_lsmfreq);
figure
nyquist(sys,sys_lsmfreq,{0.0001,pi/4});
fprintf("Open loop using LSM - frequence domain method‘s solution:\nT(time constant): %f\nL(delay): %f\n\n", 1/a1,L);

%% Open loop using LSM - frequence domain
% fprintf("Open loop using LSM - frequence domain\n")
% A=1;
% omiga = 0:pi/(len-1):pi;
% %omiga = -pi:pi/(len-1):0;
% omiga = omiga';
% len_freq = length(omiga);
% phai = [zeros(len_freq, 1) ones(len_freq, 1)];
% tao = zeros(len_freq, 1);
% 
% I_sin = zeros(len,1);
% I_cos = zeros(len,1);
% Gjw = zeros(len,1);
% 
% for i=1:len
%     for j=1:len
%         I_sin(i) = I_sin(i)+( response(j)-steady_value )*sin( omiga(i) *t(j)); % 每个wi对应的积分
%         I_cos(i) = I_cos(i)+( response(j)-steady_value )*cos( omiga(i) *t(j));
%     end
% end
% for i=1:len
% Gjw(i) = (steady_value+omiga(i)*I_sin(i)+1j*omiga(i)*I_cos(i))/A;
% end
% Gjw_mag = abs(Gjw);
% Gjw_ang = angle(Gjw);
% 
% for i=1:len
%     phai(i,1) = -omiga(i)^2*Gjw_mag(i)^2;
%     tao(i) = Gjw_mag(i)^2;
% end
% 
% theta = (inv(phai'*phai))*phai'*tao
% a1 = sqrt(theta(1))
% b0 = sqrt(theta(2))
% 
% L = ((inv(omiga'*omiga))*omiga')*(-Gjw_ang-angle(1+1j*a1*omiga))
% 
% sys_lsmfreq = tf(0.1*b0,[a1 1],'outputdelay',L) 
% response_lsmfreq = step(sys_lsmfreq,t);
% plot(t,response_lsmfreq);
% 
% fprintf("Open loop using LSM - frequence domain method‘s solution:\nT(time constant): %f\nL(delay): %f\n\n", 1/a1,L);
%% test