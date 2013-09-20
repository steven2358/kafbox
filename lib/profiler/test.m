clc
close all
clear all

N = 1000;

flops_all = zeros(N,1);
bytes_all = zeros(N,1);

kaf = krlst_profiler();
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end
    
    kaf = kaf.train_elapsed(randn,randn);
    
    flops_all(i) = kaf.lastflops();
    bytes_all(i) = kaf.lastbytes();
end
fprintf('\n');
disp(kaf)

figure; plot(bytes_all); ylabel('bytes')
figure; plot(flops_all); ylabel('flops')
% figure; semilogy(flops_all); ylabel('flops')