%% Program to post-process the data obtained from AFiD-Darcy [1]
% [1] DE PAOLI, M., Yerragolam, G. S., Lohse, D., & Verzicco, R.
%     AFiD-Darcy: A finite difference solver for numerical simulations of
%     convective porous media flows

clc ; clear ; close all

%% input parameters

% flow considered
maindirname = '../' ;   % input folder
iconfig = 1 ;   % configuration considered: 1 - RB, 2 - OS, 3 - RT
RA = 1.0e4 ;    % Rayleigh-Darcy number
Ttop = 0.0 ;    % temperature top boundary
Tbot = 1.0 ;    % temperature bottom boundary
T0 = 200 ;      % first instant in diffusive units

% results to show
i_Nusselt = 0 ;         % analyze Nusselt/flux/dissipation
i_movie_1 = 0 ;         % make animation from Results/
i_profiles_avg = 1 ;    % visualize time-averaged profiles
i_profiles_ins = 1 ;    % visualize instantaneous profiles

% parameters for plot
lw = 2 ; lw2 = 0.75 ; fs = 15 ; alpha = 1.5 ;
mygreen = [34 139 34]/255 ;
set(0,'defaultTextInterpreter','latex') ;

%% 1) Visualization of global quantities saved runtime

if i_Nusselt == 1

    % read output file
    fid = fopen(strcat(maindirname,'Results/nu_darcy.out')) ;
    data = cell2mat(textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f', ...
        'HeaderLines',1)) ;

    fclose(fid);
    time = data(:,1) ; % time
    dt = data(:,2) ;   % time step when stats are saved
    Nu_b = data(:,3) ; % dT/dx at bottom
    Nu_t = data(:,4) ; % dT/dx at top
    Nu_m = (Nu_b + Nu_t)/2 ; % mean (top-bottom) dT/dx
    chi = data(:,5) ;        % mean scalar dissipation
    c = data(:,6) ;          % mean temperature
    dcdt = data(:,7) ;       % derivative computed in AFID
    res1 = data(:,8) ;       % residual of budget first order
    c2 = data(:,9) ;         % mean squared temperature
    dc2dt = data(:,10) ;     % derivative computed in AFID
    res2 = data(:,11) ;      % residual of budget second order

    % dc2dt = gradient(c2,time) ; % d<c^2>/dt

    % for debug:
    %        time-averaged mean scalar dissipation
    % chimean_nuout = sum(dt.*chi)/(sum(dt))
    %        time-averaged Nusselt number
    % Numean_nuout = sum(dt.*(Nu_t*Ttop+Nu_b*Tbot))/(sum(dt))

    % First-order budget
    figure(101) ; clf ; box on ;
    grid off  ; hold all ; xlabel('time') ;
    set(gca,'Layer','top','fontsize',fs, ...
        'linewidth',lw2,'XScale','lin') ;
    plot(time,Nu_t,'-','linewidth',lw,'color',mygreen)
    plot(time,-Nu_b,'-','linewidth',lw,'color',[1 0 0])
    plot(time(3:end),dcdt(3:end)*RA,'-b','linewidth',lw)
    % residual of budget
    plot(time(3:end),dcdt(3:end)*RA-(-Nu_t(3:end)+Nu_b(3:end)),...
        'ok','linewidth',lw)
    plot(time,zeros(size(time)),'--y','linewidth',lw)
    hl = legend('$Nu_t$','$-Nu_b$', ...
        '$Ra\times d_t\langle C \rangle$', ...
        '$Ra\times d_t\langle C \rangle-(-Nu_t+Nu_b)$', ...
        '0') ;
    set(hl,'location','northoutside','numcol',2, ...
        'interpreter','latex','fontsize',fs)
    pbaspect([2 1 1])
    title('First-order budget')

    % Second-order budget    
    figure(102) ; clf ; box on ;
    grid off  ; hold all ; xlabel('time') ; 
    set(gca,'Layer','top','fontsize',fs, ...
        'linewidth',lw2,'XScale','lin') ;
    plot(time(3:end),-Nu_t(3:end)*Ttop,'-','linewidth',lw,'color',mygreen)
    plot(time(3:end),+Nu_b(3:end)*Tbot,'-','linewidth',lw,'color','m')
    plot(time(3:end),chi(3:end),'-b','linewidth',lw)
    plot(time(3:end),0.5*dc2dt(3:end)*RA,'-r','linewidth',lw)
    plot(time(3:end),0.5*dc2dt(3:end)*RA-...
        (-Nu_t(3:end)*Ttop+Nu_b(3:end)*Tbot-chi(3:end)),...
        '-ok','linewidth',lw)
    plot(time(3:end),zeros(size(time(3:end))),'--y','linewidth',lw)
    hl = legend('$-Nu_tT_t$','$Nu_bT_b$', ...
        '$\chi$','$Ra/2\times d_t\langle C^2 \rangle$', ...
        '$Ra/2\times d_t\langle C^2 \rangle-(-Nu_tT_t+Nu_bT_b-\chi)$','0') ;
    set(hl,'location','northoutside','numcol',2, ...
        'interpreter','latex','fontsize',fs)
    pbaspect([2 1 1])
    title('Second-order budget')

    switch iconfig

        case 2

            figure(201) ; clf ; 
            box on ; set(gca,'Layer','top','fontsize',fs,'linewidth',lw2) ;
            hold all ; grid on ;
            xlabel('$t\times Ra$') ; ylabel('$Nu_t/Ra$')
            plot(time*RA,Nu_t/RA,'o','Color',mygreen,'linewidth',lw, ...
                'MarkerFaceColor',mygreen,'Markersize',4)
            plot(time*RA,0.017*time./time,'--k','linewidth',lw)
            timesdiff = 10.^linspace(1,4,100) ;
            plot(timesdiff,1./sqrt(pi*timesdiff),'-m','linewidth',lw)
            timeshut = 10.^linspace(log10(5),6,100) ;
            Fshut = 16.8*(0.73*(timeshut-16)+31.5).^(-2) ;
            plot(timeshut*RA,Fshut,'-c','linewidth',lw)
            hl = legend('$Nu_t/Ra$','$0.017$','diffusive flux', ...
                'shutdown flux') ;
            set(hl,'location','northeast','numcol',2, ...
                'interpreter','latex','fontsize',fs)            
            box on ; pbaspect([1.52 1 1]) ;
            set(gca,"XScale",'log') ; ylim([0 0.06]) ; xlim([5e1 5e6])

        case 3

            figure(201) ; clf ; 
            box on ; set(gca,'Layer','top','fontsize',fs,'linewidth',lw2) ;
            hold all ; grid on ;
            xlabel('$t\times Ra$') ; ylabel('$\chi/Ra$')
            timesdiff = 10.^linspace(1,4,100) ;
            plot(timesdiff,1./sqrt(8*pi*timesdiff),'-k','linewidth',lw)            
            plot(time*RA,chi/RA,'ob','linewidth',lw, ...
                'MarkerFaceColor','b','Markersize',4)
            hl = legend('diffusive dissipation','simulation') ;
            set(hl,'location','northoutside','numcol',2, ...
                    'interpreter','latex','fontsize',fs)
            box on ; pbaspect([1.52 1 1]) ; 
            set(gca,"XScale",'log') ; ylim([0 0.06]) ; xlim([1e0 5e6])

    end
end


%% 2) Visualization of horizontally-averaged profiles
% (The simulation has to be concluded to generate the files)

if i_profiles_avg == 1

    % time-avergaed profiles
    dirname = strcat(maindirname,'Results/Stats/') ;
    fnameh5stats = strcat(dirname,'stafield_master.h5') ;
    info = h5info(fnameh5stats) ; % see content
    info.Datasets.Name ;
    ym = h5read(fnameh5stats,'/ym') ;
    xm = h5read(fnameh5stats,'/xm') ;
    xc = h5read(fnameh5stats,'/xc') ;
    % effective time (sum dt of each sample)
    Nsamples = h5read(fnameh5stats,'/averaging_time') ;
    % time of the simulation over which statistics are computed
    % h5read(fnameh5stats,'/averaging_interval')
    RA_read = h5read(fnameh5stats,'/Rayleigh Number') ;

    temp_avg = h5read(fnameh5stats,'/temp_m1')./Nsamples ;
    temp2_avg = h5read(fnameh5stats,'/temp_m2')./Nsamples ;
    dissth_avg = h5read(fnameh5stats,'/dissth')./Nsamples ;
    disste_avg = h5read(fnameh5stats,'/disste')./Nsamples ;
    vx_avg = h5read(fnameh5stats,'/vx_m1')./Nsamples ;
    vx2_avg = h5read(fnameh5stats,'/vx_m2')./Nsamples ;
    vy_avg = h5read(fnameh5stats,'/vy_m1')./Nsamples ;
    vy2_avg = h5read(fnameh5stats,'/vy_m2')./Nsamples ;
    vz_avg = h5read(fnameh5stats,'/vz_m1')./Nsamples ;
    vz2_avg = h5read(fnameh5stats,'/vz_m2')./Nsamples ;
    dx = xc(2:end) - xc(1:end-1) ;

    % rms of fluctuations
    temp_rms = sqrt(temp2_avg - temp_avg.^2) ;
    vx_rms = sqrt(vx2_avg - vx_avg.^2) ;
    vy_rms = sqrt(vy2_avg - vy_avg.^2) ;
    vz_rms = sqrt(vz2_avg - vz_avg.^2) ;

    figure(401) ;

    subplot(1,4,1)
    hold all ; box on ; grid on ;
    plot(temp_avg,xm,'-ok')
    xlabel('$mean(temp)$') ; ylabel('$x$')

    subplot(1,4,2)
    hold all ; box on ; grid on ;
    plot(temp_rms,xm,'-ok')
    xlabel('$rms(temp fluct)$') ; ylabel('$x$')

    subplot(1,4,3)
    hold all ; box on ; grid on ;
    plot(vx_rms,xm,'-d','Color',mygreen)
    plot(vy_rms,xm,'-sr')
    plot(vz_rms,xm,'-vb')
    xlabel('$rms(vel)$') ; ylabel('$x$')

    subplot(1,4,4)
    hold all ; box on ; grid on ;
    plot(dissth_avg,xm,'-ok')
    xlabel('$\chi$') ; ylabel('$x$')

    % show also instantaneous profiles
    if i_profiles_ins == 1

        dirname = strcat(maindirname,'Results/') ;
        fnameh5stinst = strcat(dirname,'instprof_master.h5') ;
        info = h5info(fnameh5stinst) ; % see content
        hh = info.Groups.Datasets ;
        nsamples = h5read(fnameh5stinst,'/samples');

        tot_int_time = 0 ;
        avg_int_t1 = zeros(numel(xm),1) ;
        avg_int_vx2 = zeros(numel(xm),1) ;
        avg_int_dissth = zeros(numel(xm),1) ;
        avg_int_disste = zeros(numel(xm),1) ;

        for i = 1:numel(hh)

            thiscolor = [1 0 0]*(numel(hh)-i+1)/numel(hh) ;

            deltatime = h5read(fnameh5stinst, ...
                strcat('/deltatime/',hh(i).Name));
            inst_t1 = h5read(fnameh5stinst, ...
                strcat('/temp_m1/',hh(i).Name)) ;
            inst_vx2 = h5read(fnameh5stinst, ...
                strcat('/vx_m2/',hh(i).Name)) ;
            inst_dissth = h5read(fnameh5stinst, ...
                strcat('/dissth/',hh(i).Name)) ;
            inst_disste = h5read(fnameh5stinst, ...
                strcat('/disste/',hh(i).Name)) ;

            % for debug
            avg_int_t1 = avg_int_t1 + inst_t1*deltatime ;
            avg_int_dissth = avg_int_dissth + inst_dissth*deltatime ;
            avg_int_disste = avg_int_disste + inst_disste*deltatime ;
            tot_int_time = tot_int_time + deltatime ;

            figure(601) ; sgtitle(sprintf('frame number %s',hh(i).Name))

            subplot(1,2,1) ; hold all
            plot(inst_t1,xm,'-','Color',thiscolor) ;
            box on ; grid on ; set(gca,'layer','top')
            ylim([min(xm) max(xm)]) ;
            xlabel('temp') ; ylabel('$x$')

            subplot(1,2,2) ; hold all
            plot(inst_vx2,xm,'-','Color',thiscolor) ;
            box on ; grid on ; set(gca,'layer','top')
            ylim([min(xm) max(xm)]) ;
            xlabel('$V_x^2$') ; ylabel('$x$')

            drawnow

        end

    end

end

%% 3) Visualization of intantaneous scalar fields

if i_movie_1 == 1

    % read files
    outdir = strcat(maindirname,'movies/') ;
    if ~exist(outdir, 'dir') ;  mkdir(outdir) ; end
    dirname = strcat(maindirname,'Results/Temp/') ;
    % Replace "Temp" with "Vx", "Vy", "Vz" to visualize the other fields

    % read and prepare grid
    fnameh5grid = strcat(maindirname,'Results/','cordin_info.h5') ;
    %xm = h5read(fnameh5grid,'/xm') ;
    ym = h5read(fnameh5grid,'/ym') ;
    %zm = h5read(fnameh5grid,'/zm') ;
    xc = h5read(fnameh5grid,'/xc') ;
    %yc = h5read(fnameh5grid,'/yc') ;
    %zc = h5read(fnameh5grid,'/zc') ;
    [X,Y] = meshgrid(ym,xc) ;

    filedet = dir(dirname) ;
    vidfile = VideoWriter(strcat(outdir,'temperature.mp4'),'MPEG-4');
    vidfile.FrameRate = 5; % fps
    open(vidfile);

    % start from 3 to include first frame

    fprintf('Processing movie  ... \n ')
    for i = 3:numel(filedet)

        fnameh5 = strcat(dirname,filedet(i).name) ;
        if strcmp(fnameh5,strcat(dirname,'.DS_Store'))==1
            delete(fnameh5) ; % delete macOS hidden files
        else
            % hinf = h5info(fnameh5); % get file info
            % disp(hinf) ;
            temp = h5read(fnameh5,'/var') ;

            f = figure('visible','off') ; clf ; box on ; colormap(hot(512))
            surf(X,Y,temp) ; shading('interp'); view(2) ; c=colorbar ;
            set(gca,'Layer','top') ; grid off  ; box on
            h = get(gca,'DataAspectRatio') ;
            if h(3)==1
                set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
            else
                set(gca,'DataAspectRatio',[1 1 h(3)])
            end
            set(gca,'FontSize',round(15/max(ym)), ...
                'LineWidth',ceil(2/max(ym)))
            set(c,'Fontsize',round(15/max(ym)), ...
                'TickLabelInterpreter','latex', ...
                'LineWidth',ceil(2/max(ym)))
            % clim([0 1]) ; zlim([0 1])
            % small portion of the domain is visualized
            if iconfig == 2 ;  xlim([0 2]) ; end
            if iconfig == 3 ;  xlim([0 1]) ; end
            frame = getframe(gcf) ;
            writeVideo(vidfile,frame) ;
        end

    end

    close(vidfile)

end



