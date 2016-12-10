function [K] = BPF_2D(x) 
F = zeros(size(x,1),size(x,2));
for i = 1:size(x,1)
    F(i,:) = BPF_1D(x(i,:)');
for i = 1:size(x,2)
    K(:,i) = BPF_1D(F(:,i))'; 
end
end

