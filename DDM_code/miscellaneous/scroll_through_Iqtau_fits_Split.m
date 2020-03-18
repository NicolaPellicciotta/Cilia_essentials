function [] = scroll_through_Iqtau_fits_Split( Iqtau, fits_array, tau_array )
%Iqtau_reader simple GUI to scroll through the modes of an Iqtau matrix
%   Detailed explanation goes here

hf = figure;
hf.Position = hf.Position .*[1 1 2 1];

%% exit button
STOP = 0;
STOP_ = uicontrol('style','togglebutton','string','CLOSE','min',0,'max',1,'value',0);
set(STOP_,'units','normalized','position',[.45 .0 .1 .05]);

%% mode slider
mode_ = uicontrol('style','slider','min',1,'max',size(Iqtau,1),'value',1);
set(mode_,'units','normalized','position',[.10 .10 .80 .03]);
uicontrol('style','text','string','mode','units','normalized','position',[.03 .10 .06 .03]);

%% initial image

ha = gca;
hold on;
ha.Position = [.05 .25 .90 .70];
%  ha.XScale = 'log';
%  ha.YScale = 'log';
mode = get(mode_,'value');
hp1 = plot(tau_array(),Iqtau(mode,:),'b.');  %-fits_array{57}(tau_array(end)
hp2 = plot(tau_array,fits_array{mode}(tau_array),'r');
    hp3 = plot(tau_array([1:7:end]),Iqtau(mode,[1:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp4 = plot(tau_array([2:7:end]),Iqtau(mode,[2:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp5 = plot(tau_array([3:7:end]),Iqtau(mode,[3:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp6 = plot(tau_array([4:7:end]),Iqtau(mode,[4:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp7 = plot(tau_array([5:7:end]),Iqtau(mode,[5:7:end]),'g-');  %-fits_array{57}(tau_array(end)
        hp8 = plot(tau_array([6:7:end]),Iqtau(mode,[6:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp9 = plot(tau_array([7:7:end]),Iqtau(mode,[7:7:end]),'-');  %-fits_array{57}(tau_array(end)


xlabel('Lag time \tau / s');
hyl = ylabel('Squared Amplitude of the Fourier mode/a.u.');
hyl.Units = 'Normalized';
hyl.Position = hyl.Position - [0.01 0 0];

%% adjust and show
while STOP==0
    
    new_mode = round(get(mode_,'value'));
    
    if(new_mode ~= mode)
        mode = new_mode;
        uicontrol('style','text','string',num2str(mode),'units','normalized','position',[.91 .10 .06 .03]);
        delete(hp1);
        delete(hp2);
        delete(hp3);
        delete(hp4);
        delete(hp5);
        delete(hp6);
        delete(hp7);
        delete(hp8);
        delete(hp9);
        
           
%         hp = plot(Iqtau(mode,1:floor(end/2)),'b');
        hp1 = plot(tau_array(),Iqtau(mode,:),'b.');
        hold on;
        hp2 = plot(tau_array,fits_array{mode}(tau_array),'r');
        hp3 = plot(tau_array([1:7:end]),Iqtau(mode,[1:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp4 = plot(tau_array([2:7:end]),Iqtau(mode,[2:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp5 = plot(tau_array([3:7:end]),Iqtau(mode,[3:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp6 = plot(tau_array([4:7:end]),Iqtau(mode,[4:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp7 = plot(tau_array([5:7:end]),Iqtau(mode,[5:7:end]),'g-');  %-fits_array{57}(tau_array(end)
        hp8 = plot(tau_array([6:7:end]),Iqtau(mode,[6:7:end]),'-');  %-fits_array{57}(tau_array(end)
        hp9 = plot(tau_array([7:7:end]),Iqtau(mode,[7:7:end]),'-');  %-fits_array{57}(tau_array(end)
%         hyl.Position = hylPosition;
    end
    pause(.05);
    
    STOP = get(STOP_,'value');
end

close(hf)

end

