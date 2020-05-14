%function analysis_spec()
%{
    JASON CARDARELLI
    CODE TO ANALYZE AND VISUALIZE OUTPUT FROM MAG SPEC CODE
    DATA FORMAT:
        ROWS: PARTICLES
        COLS: TIME STEPS
%}
close all
    posx = readmatrix("XPOS.csv",'NumHeaderLines',0);
    posy = readmatrix("YPOS.csv",'NumHeaderLines',0);
    posz = readmatrix("ZPOS.csv",'NumHeaderLines',0);
    velx = readmatrix("MOMENTUM_X.csv",'NumHeaderLines',0);
    vely = readmatrix("MOMENTUM_Y.csv",'NumHeaderLines',0);
    velz = readmatrix("MOMENTUM_Z.csv",'NumHeaderLines',0);
    time = readmatrix("TIME.csv",'NumHeaderLines',0);
    energy = readmatrix("ENERGY.csv",'NumHeaderLines',0);
    magnet = readmatrix("MAGNETS.csv"); %index, Bz, x, y, z, length, width
        mag_col_names = {'Magnet_Index','Bz','pos_x','pos_y','pos_z','length','width','height'};
        magnet = array2table(magnet,'VariableNames',mag_col_names);
        clear mag_col_names;
    screen = readmatrix("SCREENS.csv");
        screen_col_names = {'Screen_index','edgeX','edgeY','edgeZ','angle','length','height'};
        screen = array2table(screen, 'VariableNames', screen_col_names);
        clear screen_col_names;
    
        
    [num_par,~]    = size(posx);
    [num_mag,~]    = size(magnet);
    [num_screen,~] = size(screen);
    
    
    normalizing_omega = 1.602*10^(-19)*1/(1000*9.11*10^(-31)/0.511);
    c = 3*10^8;
    central_en_x = 1000;
    
    
        %initial conditions:
    int_vel = [velx(:,1),vely(:,1),velz(:,1)]; 
    int_pos = [posx(:,1),posy(:,1),posz(:,1)];
        
        %final coniditons for each particle: 
    posx_nonan = ~isnan(posx);
    posx_Indices = arrayfun(@(x) find(posx_nonan(x, :), 1, 'last'), 1:num_par);
    posx_values = arrayfun(@(x,y) posx(y,x), posx_Indices, 1:num_par);
    posy_nonan = ~isnan(posy);
    posy_Indices = arrayfun(@(x) find(posy_nonan(x, :), 1, 'last'), 1:num_par);
    posy_values = arrayfun(@(x,y) posy(y,x), posy_Indices, 1:num_par);
    posz_nonan = ~isnan(posz);
    posz_Indices = arrayfun(@(x) find(posz_nonan(x, :), 1, 'last'), 1:num_par);
    posz_values = arrayfun(@(x,y) posz(y,x), posz_Indices, 1:num_par);
    fin_pos = [posx_values.',posy_values.',posz_values.'];
    clear posx_nonan posy_nonan posz_nonan 
    clear posx_Indices posy_Indices posz_Indices
    clear posx_values posy_values posz_values
    velx_nonan = ~isnan(velx);
    velx_Indices = arrayfun(@(x) find(velx_nonan(x, :), 1, 'last'), 1:num_par);
    velx_values = arrayfun(@(x,y) velx(y,x), velx_Indices, 1:num_par);
    vely_nonan = ~isnan(vely);
    vely_Indices = arrayfun(@(x) find(vely_nonan(x, :), 1, 'last'), 1:num_par);
    vely_values = arrayfun(@(x,y) vely(y,x), vely_Indices, 1:num_par);
    velz_nonan = ~isnan(velz);
    velz_Indices = arrayfun(@(x) find(velz_nonan(x, :), 1, 'last'), 1:num_par);
    velz_values = arrayfun(@(x,y) velz(y,x), velz_Indices, 1:num_par);
    fin_vel = [velx_values.',vely_values.',velz_values.'];
    clear velx_nonan vely_nonan velz_nonan 
    clear velx_Indices vely_Indices velz_Indices
    clear velx_values vely_values velz_values   
        
   
    
    figure
    hold on
    for ii = 1:num_mag
        %rectangle('Position',[magnet.pos_x(ii).*(c/normalizing_omega) (magnet.pos_y(ii)-(magnet.width(ii)/2)).*(c/normalizing_omega) magnet.length(ii).*(c/normalizing_omega) magnet.width(ii).*(c/normalizing_omega)])
        rectangle('Position',[magnet.pos_x(ii) (magnet.pos_y(ii)-(magnet.width(ii)/2)) magnet.length(ii) magnet.width(ii)])
    end
    
    for ii = 1:num_screen
        screen_x = [screen.edgeX(ii), screen.edgeX(ii)+(screen.length(ii)*cosd(screen.angle(ii)))];
        %screen_x = screen_x.*(c/normalizing_omega);
        screen_y = [screen.edgeY(ii), screen.edgeY(ii)+(screen.length(ii)*sind(screen.angle(ii)))];
        %screen_y = screen_y.*(c/normalizing_omega);
        line(screen_x,screen_y);
    end
    
    for ii = 1:num_par
        plot(posx(ii,:),posy(ii,:),'DisplayName',['Initial KE_X = ', num2str(energy(ii,1).*0.511), ' MeV'])  
        %plot(posx(ii,:).*(c/normalizing_omega),posy(ii,:).*(c/normalizing_omega),'DisplayName',['Initial KE_X = ', num2str(energy(ii,1)), ' MeV'])
    end
    axis equal
    legend
    %xlabel('Longitudinal Position [x \omega_c/c]')
    xlabel('Longitudinal Position [m]')
    ylabel('Laterial Position [m]')
    %title([num2str(num_par),' Particle Trajectories Through ', num2str(num_mag), ' Magnets of Strength ', num2str(table2array(magnet(1,2))), ' and ', num2str(table2array(magnet(2,2))), ' q_eB/m_e\omega_c'])
    hold off
    
    %{
    figure(2)
    hold on
    
    hold off
    %}
%end
