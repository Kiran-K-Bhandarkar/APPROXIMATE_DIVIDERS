clc;
clear;
fb = zeros(1, 256);
Xb = zeros(1, 256);
for i= 0:255
    i_bin = dec2bin(i, 8);
    decimal_val = (str2double(i_bin(1))/2) + (str2double(i_bin(2))/4) + ...
    (str2double(i_bin(3))/8) + (str2double(i_bin(4))/16) + ...
    (str2double(i_bin(5))/32) + (str2double(i_bin(6))/64) + ...
    (str2double(i_bin(7))/128) + (str2double(i_bin(8))/256);
    fb(1,i+1)=decimal_val;
end

for i=1:length(fb)
    Xb(i) = 1/(1+fb(i));
end

figure;
subplot(2,2,1);
axis([0 1 0.4 1]);
title('Xb vs fb (Exact)'),
xlabel('fb');
ylabel('Xb');
hold on;
grid on;
plot(fb, Xb);

% Two slices
Xb_d = [Xb(1) Xb(128) Xb(256)];
fb_d = [fb(1) fb(128) fb(256)];
% % fb = 0 to 0.4961
% slope = -0.2345;
% intercept = 0.8342;

% % fb = 0.4961 to 0.9961
% slope =  -0.1184;
% intercept =  0.5847;

subplot(2,2,2);
axis([0 1 0.4 1]);
title('Xb vs fb (r = 2)'),
xlabel('fb');
ylabel('Xb');
hold on;
grid on;
plot(fb, Xb, fb_d,Xb_d);

% Four slices
Xb_d = [Xb(1) Xb(64) Xb(128) Xb(192) Xb(256)];
fb_d = [fb(1) fb(64) fb(128) fb(192) fb(256)];
subplot(2,2,3);
axis([0 1 0.4 1]);
title('Xb vs fb (r = 4)'),
xlabel('fb');
ylabel('Xb');
hold on;
grid on;
plot(fb, Xb, fb_d,Xb_d);

% Eight slices
Xb_d = [Xb(1) Xb(32) Xb(64) Xb(96) Xb(128) Xb(160) Xb(192) ...
    Xb(224) Xb(256)];
fb_d = [fb(1) fb(32) fb(64) fb(96) fb(128) fb(160) fb(192) ...
    fb(224) fb(256)];
subplot(2,2,4);
axis([0 1 0.4 1]);
title('Xb vs fb (r = 8)'),
xlabel('fb');
ylabel('Xb');
hold on;
grid on;
plot(fb, Xb, fb_d,Xb_d);