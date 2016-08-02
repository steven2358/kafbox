% Demonstration of the basic operation of a kernel adaptive filtering
% algorithm with fixed budget. This script generates several images and
% saves them to disk.
%
% KRLS-T is chosen as the example algorithm. Inspired by
% http://www.tsc.uc3m.es/~miguel/MLG/adjuntos/slidesKRLST.pdf
%
% This file is part of the Kernel Adaptive Filtering Toolbox for Matlab.
% https://github.com/steven2358/kafbox/

close all
clear
rng('default'); rng(11);

%% PARAMETERS

n = 10; % number of training data
n_star = 201; % number of points to evaluate

ell = 2; % length scale
sigvar = 1; % signal variance
noisevar = .02; % noise variance

budget = 5; % KRLS-T budget

savefig = 1;
manual = 0;

%% GENERATE DATA

fprintf('Generating data...\n');
covfunc = {'kafbox_covSum', {'kafbox_covSEiso','kafbox_covNoise'}};
loghyper = [log(ell); log(sigvar); log(noisevar)];

x = (1:n)' + 0.2*randn(n,1);
y = chol(feval(covfunc{:}, loghyper, x))'*randn(n,1);      % Cholesky decomp.

x_star = linspace(0,n+2,n_star)';

covfunc_no_noise = {'kafbox_covSEiso'};
loghyper_no_noise = [log(ell); log(sigvar)];

%% RUN ALGORITHM

kaf = krlst_split(struct('kerneltype','gauss','kernelpar',ell,...
    'sn2',noisevar/sigvar,'M',budget+1));

fprintf('Running KRLS-T\n');
t1 = tic;
for i=1:n,
    if ~mod(i,floor(n/10)), fprintf('.'); end
    Qold = kaf.Q;
    kaf = kaf.train(x(i),y(i));
    
    for stage = 1:2
        switch stage
            case 1 % stage 1: after including a sample
                stage_string = 'after inclusion';
                fname_i = sprintf('fig/sample_run_%02da',i);
            case 2 % stage 2: after pruning
                if length(kaf.mu) > budget
                    kaf = kaf.prune(Qold,1);
                    stage_string = 'after pruning  ';
                    fname_i = sprintf('fig/sample_run_%02db',i);
                end
        end
        
        y_est = kaf.evaluate(x(1:i));
        [y_star_est,y_star_var] = kaf.evaluate(x_star);
        
        y_star_var = y_star_var/kaf.s02*sigvar; % fix scale
        
        x_dict = kaf.dict;
        y_dict = kaf.evaluate(x_dict);
        
        papersize = [8 4];
        
        figure(1); clf
        set(gcf,'Position',[500 400 papersize(1)*100 papersize(2)*100])
        hold all
        f = [y_star_est+2*sqrt(y_star_var); flip(y_star_est-2*sqrt(y_star_var),1)];
        p0 = fill([x_star; flip(x_star,1)], f, [7 7 7]/8,'EdgeColor','none');
        p1 = plot(x_star,y_star_est,'r','LineWidth',1.5);
        p2 = plot(x_dict,y_dict,'ro','LineWidth',2,'MarkerSize',8);
        p3 = plot(x(1:i),y(1:i),'bx','LineWidth',2,'MarkerSize',8);
        
        set(gca,'box','on');
        set(gca, 'FontSize', 14);
        set(gca, 'FontName', 'Helvetica');
        
        axis([0 n+2 floor(min(y))-1 ceil(max(y))+1])
        legend([p3,p1,p2],'data','prediction','dictionary bases')
        
        title(sprintf('KRLS-T iteration %d: %s',i,stage_string));
        
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', papersize);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 papersize(1) papersize(2)]);
        
        style = hgexport('factorystyle');
        style.Bounds = 'tight';
        hgexport(figure(1),'-clipboard',style,'applystyle', true);
        set(gcf, 'renderer', 'painters');
        
        drawnow
        if manual
            keyboard; %#ok<UNRCH>
        end
        
        if savefig
            if ~exist('fig','file')
                mkdir('fig');
            end
            % write individual png
            print('-dpng','-r100',fname_i);
            % write animated gif
            fname = 'fig/sample_run.gif';
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i==1 && stage==1
                imwrite(imind,cm,fname,'gif','Loopcount',inf,'DelayTime',1.5);
            else
                if stage==2 && i==n && kaf.pruned
                    imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',2.5);
                elseif ~(stage==2 && ~kaf.pruned)
                    imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',1.5);
                end
            end
        end
    end
    
end
fprintf('\n')
