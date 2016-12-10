function [F] = BPF_1D(x)
a = 5;
N = size(x,1)-1;% N отсчетов на апертуре
x1 = x(1:size(x)/2);
x2 = x(fix(size(x)/2+1):size(x)-1);
x1 = cat(1,zeros(500,1),x1);
x2 = cat(1,x2,zeros(500,1));
f = cat(1,x2,x1);
F = fft(f);
%F = fftshift(F); 
F = F*2*a/N;
f1 = F(1:size(F)/2);
f2 = F(size(F)/2+1:size(F));
F = cat(1,f2,f1);
F = F(size(F)/2-N/2:size(F)/2+N/2);
end

