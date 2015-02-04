classdef Visualize < handle
    properties
        hf1;
        hf2;
        hf3;
    end
    methods
        function Initialize(obj)
            obj.hf1 = figure;
            set(obj.hf1, 'position', [50, 50, 640, 640]);
            obj.hf2 = figure;
            set(obj.hf2, 'position', [200, 80, 800, 600]);
            obj.hf3 = figure;            
            set(obj.hf3, 'position', [220, 100, 800, 600]);
        end
        function PlotVelocity(obj, SimulSys)
            pb = SimulSys.pb;
%             clf(obj.hf1);
%             clf(obj.hf2);
%             clf(obj.hf3);
            %% Plot Particles
            figure(obj.hf1);
            clf
            hold on
            plot_particles(obj.hf1, SimulSys.part, SimulSys.ghost, SimulSys.pb);
            %%
            pb = SimulSys.pb;
            loc = [2*pb.h; -2.5*pb.h];
            plot_neighbours(obj.hf1, loc, pb, SimulSys.part, SimulSys.ghost);
            plot(loc(1), loc(2), 'r*', 'markerfacecolor', 'r');
            ip = round(pb.nx/2);
            jp = pb.ny-1;
            loc = SimulSys.part.r(:,(ip-1)*pb.ny + jp);
            plot_neighbours(obj.hf1, loc, pb, SimulSys.part, SimulSys.ghost);
            plot(loc(1), loc(2), 'r*', 'markerfacecolor', 'r');
            hold off
            %% Plot Pressure
            figure(obj.hf2);
            x_lim = [0 pb.L];
            y_lim = [-pb.b pb.b];

            x_min = x_lim(1) - 2*pb.h;
            x_max = x_lim(2) + 2*pb.h;

            y_min = y_lim(1) - 2*pb.h;
            y_max = y_lim(2) + 2*pb.h;

            subplot(3,1,1)
            x = SimulSys.part.r(1,:)';
            y = SimulSys.part.r(2,:)';
            z = SimulSys.part.p';
            tri=delaunay(x,y);
            h = trisurf(tri, x, y, z);
            l = light('Position',[-50 -15 29]);
            lighting phong
            shading interp
            view(0, 90);
            xlim([x_min x_max]);
            ylim([y_min y_max]);
            colorbar

            subplot(3,1,2)
            x=[SimulSys.ghost.r(1,:)' ; SimulSys.part.r(1,:)'];
            y=[SimulSys.ghost.r(2,:)' ; SimulSys.part.r(2,:)'];
            z=[SimulSys.ghost.p' ; nan*ones(pb.N,1)];
            tri=delaunay(x,y);
            h = trisurf(tri, x, y, z);
            %axis vis3d
            l = light('Position',[-50 -15 29]);
            lighting phong
            shading interp
            view(0, 90)
            xlim([x_min x_max]);
            ylim([y_min y_max]);
            colorbar

            subplot(3,1,3)
            x=[SimulSys.ghost.r(1,:)' ; SimulSys.part.r(1,:)'];
            y=[SimulSys.ghost.r(2,:)' ; SimulSys.part.r(2,:)'];
            z=[SimulSys.ghost.p' ; SimulSys.part.p'];
            tri=delaunay(x,y);
            h = trisurf(tri, x, y, z);
            %axis vis3d
            l = light('Position',[-50 -15 29]);
            lighting phong
            shading interp
            view(0, 90)
            xlim([x_min x_max]);
            ylim([y_min y_max]);
            colorbar
            %% Plot Velocity
            figure(obj.hf3);
            subplot(2,1,1)
            x=[SimulSys.ghost.r(1,:)' ; SimulSys.part.r(1,:)'];
            y=[SimulSys.ghost.r(2,:)' ; SimulSys.part.r(2,:)'];
            z=[SimulSys.ghost.v(1,:)' ; SimulSys.part.v(1,:)'];
            tri=delaunay(x,y);
            h = trisurf(tri, x, y, z);
            %axis vis3d
            l = light('Position',[-50 -15 29]);
            lighting phong
            shading interp
            view(0, 90)
            xlim([x_min x_max]);
            ylim([y_min y_max]);
            colorbar

            subplot(2,1,2)
            x=[SimulSys.ghost.r(1,:)' ; SimulSys.part.r(1,:)'];
            y=[SimulSys.ghost.r(2,:)' ; SimulSys.part.r(2,:)'];
            z=[SimulSys.ghost.v(2,:)' ; SimulSys.part.v(2,:)'];
            tri=delaunay(x,y);
            h = trisurf(tri, x, y, z);
            %axis vis3d
            l = light('Position',[-50 -15 29]);
            lighting phong
            shading interp
            view(0, 90)
            xlim([x_min x_max]);
            ylim([y_min y_max]);   
            colorbar
        end
    end % methods
end % class
            
            