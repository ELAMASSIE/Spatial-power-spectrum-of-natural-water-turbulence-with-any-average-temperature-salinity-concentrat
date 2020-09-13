function Ot = Y_Schmidt(T,uT,S,uS)
%salinity Schmidt number varies with av. temperature and salinity
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
    if ~isequal((T<273.15)+(T>453.15),zeros(size(T)))
        warning('Temperature is out of range for Viscosity function 0<T<180 C');
    end

    if ~isequal((S<0)+(S>150),zeros(size(S)))
        warning('Salinity is out of range for Viscosity function 0<S<150 g/kg');
    end

    %% BEGIN
    Ot = SW_Kviscosity(T,'K',S,uS)./Y_Sdif(T,'K',S,uS) * 10^10;

end

