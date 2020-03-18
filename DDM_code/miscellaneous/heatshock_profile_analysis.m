close all

filename = 'D:\Data\Cilia\Data\2015_03_04\killing_cilia_profile.04Mar2015_17.17.51.movie';
px2mum = 0.29; %um/px Grasshopper 20X

N_frames_to_std = 20; %number of frames on which to do std



mo = moviereader(filename);

T_tot = floor(mo.NumberOfFrames/mo.FrameRate); %total time, given by number of frames and framerate. in seconds

time_vec = 0:0.15:T_tot-1; %vector of time, from 0 to 1s before the video ends.

first_frame = mo.read(1); %select roi on first frame
imshow(first_frame,[]);

rect_out = getrect(gca);    %nifty function to convert rectangle to rows and columns

[rows,cols] = rect2sub(rect_out);

imshow(first_frame(rows,cols),[]);
drawnow;

%%
%for every entry of time_vec load some frames, do std
std_stack = zeros(length(rows),length(cols),length(time_vec));
fprintf('\n');
for i = 1:length(time_vec);
    
    t = time_vec(i);
    Fi = round(t*mo.FrameRate +1); %initial frame to load
    Ff = Fi + N_frames_to_std;  %final frame to load
    fs = mo.read([Fi, Ff]);
    fs = fs(rows,cols,:);
    std_stack(:,:,i) = std(single(fs),1,3);
        
    if mod(i,10) == 0
        fprintf('\b\b\b\b\b\b\b\b\b%.4d/%.4d',i,length(time_vec))
    end
    
end
%%
figure
std_kimograph = squeeze(sum(std_stack,1)); %fake kimograph (actually sum of each column)
imshow(imadjust(mat2gray(std_kimograph')),[]);
ha = gca;
ha.Visible = 'on';
% ha.Box = 'off';
ha.YTickLabel = time_vec(ha.YTick);
ha.TickLength = ha.TickLength./2;
hyl = ylabel('Time, [s]');
hxl = xlabel('Pixels');



%%
[x, t] = ginputc('Color','r','LineWidth',0.2, 'LineStyle', ':', 'ShowPoints', true, 'ConnectPoints', false);
t = round(t);
%%
fit_out = fit(x,time_vec(t)',fittype('a*(x-x0)^2 + t0','independent','x', 'options', fitoptions('Method','NonlinearLeastSquares', 'StartPoint', [mean(x) 0  min(time_vec(t))] )));

xplot = linspace(x(1),x(end), length(x)*100);

figure
plot(x, time_vec(t)', 'o');
hold on
plot(xplot, fit_out(xplot));
xlabel('Pixels');
ylabel('Time, [s]');
legend('data','parabolic fit','Location','NorthEast');


%%
xx = abs(x-fit_out.x0);
time = time_vec(t) - fit_out.t0; %now tt=a*xx^2 => xx = sqrt(tt/a)
xxplot = abs(xplot - fit_out.x0);

figure
plot(time',xx*px2mum, 'o');
hold on
plot(fit_out.a.*xxplot.^2,xxplot*px2mum);
ylabel('Radius of ciliary inactivity, [\mum]');
xlabel('Time laser on, [s]');
hl = legend('data','fit \propto time^{1/2}');
hl.Location = 'SouthEast';