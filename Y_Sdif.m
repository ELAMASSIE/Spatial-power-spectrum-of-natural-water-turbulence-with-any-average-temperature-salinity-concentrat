function Ot = Y_Sdif(T,uT,S,uS)
%salinity diffusion coefficient [E-10 m^2/s)] varies with av. temperature and salinity
 %% CHECK INPUT ARGUMENTS

    % CHECK THAT S&T HAVE SAME SHAPE
    if ~isequal(size(S),size(T))
        error('check_stp: S & T must have same dimensions');
    end

    % CONVERT TEMPERATURE INPUT TO °C
    switch lower(uT)
        case 'c'
            T = T + 273.15;
        case 'k'
            T = T;
        otherwise
            error('Not a recognized temperature unit. Please use ''C'', ''K''');
    end

    % CONVERT SALINITY TO PPT
    switch lower(uS)
        case 'ppt'
        case 'ppm'
            S = S./1000;
        case 'w'
            S = S.*1000;
        case '%'
            S = S.*10;
        otherwise
            error('Not a recognized salinity unit. Please use ''ppt'', ''ppm'', ''w'', or ''%''');
    end

    % CHECK THAT S & T ARE WITHIN THE FUNCTION RANGE
    if ~isequal((T<273.15)+(T>273.15+180),zeros(size(T)))
        warning('Temperature is out of range for Viscosity function 0<T<180 C');
    end

    if ~isequal((S<0)+(S>150),zeros(size(S)))
        warning('Salinity is out of range for Viscosity function 0<S<150 g/kg');
    end

    %% BEGIN
    data_T = [15 20 25 27 30 25 25 25 25 25]'+ 273.15;
    data_S = [34.9 34.9 34.9 34.9 34.9 0 8.9 17.7 26.3 34.9]';
    data_d = [12.86 14.5 17.71 17.37 18.46 20.3	18.71 18.38 17.96 17.71]';
    
    data_v = SW_Viscosity(data_T,'K',data_S,'ppt'); 
    % Direct average
%     data_M = data_d./data_T.*data_v;
%     Ot = T./SW_my_Viscosity(T,'K',S,'ppt').*mean(data_M);

    % least square method  data_T./data_v * k = data_d
    k = (mean(data_T./data_v.*data_d) - mean(data_T./data_v)*mean(data_d))/(mean((data_T./data_v).^2) - (mean(data_T./data_v)).^2);
    Ot = T./SW_Viscosity(T,'K',S,'ppt').*k;

end

