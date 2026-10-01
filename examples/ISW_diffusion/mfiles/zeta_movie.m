addpath '/Users/brookepauken/suntans/suntans-master/mfiles';
close all;

ANIMATE=~false;
PLOTZETA=1;
PLOTOMEGA=0;
PLOTOMEGAALT=1;
km = 1000;

modifier='new_adv_test_QUICK';
dirname = '../data_oscitest_newadv';

% dirname='~/oscitests/copy_hybrid_10T_800_uniform_qs_layerdiff800';
% modifier='hybrid_10T_800_uniform_qs_layerdiff800';

% modifier='hybrid_10T_800tests_dzf6_qs_05dt800';
% dirname = '/Users/brookepauken/oscitests/copy_hybrid_10T_800tests_dzf6_qs_05dt800';

suntansfile = [dirname,'/suntans.dat'];

c = load([dirname,'/cells.dat'],'-ascii');

xv = c(:,2);
Nc = length(xv);


Nkmax = getvalue(suntansfile,'Nkmax');
nsteps = getvalue(suntansfile,'nsteps');
ntout = getvalue(suntansfile,'ntout');
dt = getvalue(suntansfile,'dt');

[XV, ~] = meshgrid(xv, 1:(Nkmax+1));
XV=XV';

data = fread(fopen([dirname,'/vert_dzz.dat'],'rb'),'float64');
nout = length(data)/Nkmax/Nc;

if(PLOTZETA)
    vidtit=['/Users/brookepauken/suntans-gvc-mild-slope-2025-results/quick-isw-test-results/zetaplot_', modifier];
else
    vidtit=['/Users/brookepauken/suntans-gvc-mild-slope-2025-results/quick-isw-test-results/omegzetaplotlog_', modifier];
end
vid=VideoWriter(vidtit, 'MPEG-4');
vid.FrameRate = 30;
%vid.FrameRate = 20;
vid.Quality = 90;
open(vid)

for n=1:nout
    fprintf('n = %d of %d\n',n,nout);
    zeta = get_zeta(dirname,n);
    

    if(ANIMATE)
        if(PLOTZETA)
            figure(1);  
            cla;
            plot(xv/km,zeta,'k-');
            axis tight;
            ylim([-305, 0])
            drawnow;
        elseif(PLOTOMEGA)
            figure(1);  
            cla;
            omega = get_omega(dirname,n);
            pcolor(XV./km, zeta, omega);
            shading flat
            axis tight;
            c = colorbar;
            c.Label.String = 'omega (m/s)';
            c.Label.Interpreter = 'latex';
            %clim([-6e-5, 5e-5]);
            hold on
            plot(xv/km,zeta,'w-');
            xlabel('x (km)')
            ylabel('z (m)')
            
            %ylim([-305, 0])
            drawnow;
        elseif(PLOTOMEGAALT)
            figure(1);  
            cla;
            omega = get_omega(dirname,n);

            % make a custom colormap
            % n = 256; % how many levels do you want?
            % cvu = linspace(0,1,n).';
            % cvd = flipud(cvu);
            % cmap = [cvu cvu cvd];

            % instead of using plot(), use surf()
            % because it can do color interpolation
            for(i=1:Nkmax+1)
                x = xv;
                y = zeta(:, i);
                z = log10(abs(omega(:, i)));
                h = surf([x(:) x(:)],[y(:) y(:)],[z(:) z(:)]);

            
                %drawnow
                set(h,'facecolor','none','edgecolor','interp');
                view(2);
                %set(h,'linewidth',3);
                %pcolor(XV./km, zeta, omega);
                %shading flat
                % axis tight;

                clim([-10, -3]);
                hold on
            end

            colormap('parula')
            c = colorbar;
            c.Label.String = 'log abs omega';
            c.Label.Interpreter = 'latex';
            grid off
            %pause;
            % plot(xv/km,zeta,'w-');
            % xlabel('x (km)')
            % ylabel('z (m)')
            
            ylim([-300, 0])
            drawnow;

        end

        F = getframe(figure(1)); 
        writeVideo(vid, F);
    end
end

% figure;
% [XV, ~] = meshgrid(xv, 1:(Nkmax+1));
% surf(XV', zeta, omega)

close(vid);

