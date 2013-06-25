clc
clear all

N = 1000;

flops_all = zeros(N,1);
bytes_all = zeros(N,1);

kaf = rls_profiler(struct('lambda',.999,'c',1E-4));
for i=1:N,
    if ~mod(i,floor(N/10)), fprintf('.'); end
    
    kaf = kaf.train_elapsed(rand,rand);
    
    flops_all(i) = kaf.lastflops();
    bytes_all(i) = kaf.lastbytes();
end
disp(kaf)

figure;plot(bytes_all)
figure;plot(flops_all)