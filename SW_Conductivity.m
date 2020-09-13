function k = SW_Conductivity(T,uT,S,uS)
    % SW_Conductivity    Thermal conductivity of seawater
    %=========================================================================
    % USAGE:  k = SW_Conductivity(T,uT,S,uS)
    %
    % DESCRIPTION:
    %   Thermal conductivity of seawater at 0.1 MPa given by [1]
    %   Values at temperature higher than the normal boiling temperature
    %   are calculated at the saturation pressure.
    %
    % INPUT:
    %   T  = temperature
    %   uT = temperature unit
    %        'C'  : [degree Celsius] (ITS-90)
    %        'K'  : [Kelvin]
    %        'F'  : [degree Fahrenheit]
    %        'R'  : [Rankine]
    %   S  = salinity
    %   uS = salinity unit
    %        'ppt': [g/kg] (reference-composition salinity)
    %        'ppm': [mg/kg] (in parts per million)
    %        'w'  : [kg/kg] (mass fraction)
    %        '%'  : [kg/kg] (in parts per hundred)
    %
    %   Note: T and S must have the same dimensions
    %
    % OUTPUT:
    %   k = thermal conductivity [W/m K]
    %
    %   Note: k will have the same dimensions as T and S
    %
    % VALIDITY: 0 < T < 180 C; 0 < S < 160 g/kg
    %
    % ACCURACY: 3.0%
    %
    % REVISION HISTORY:
    %   2009-12-18: Mostafa H. Sharqawy (mhamed@mit.edu), MIT
    %               - Initial version
    %   2012-06-06: Karan H. Mistry (mistry@alum.mit.edu), MIT
    %               - Allow T,S input in various units
    %               - Allow T,S to be matrices of any size
    %
    % DISCLAIMER:
    %   This software is provided "as is" without warranty of any kind.
    %   See the file sw_copy.m for conditions of use and licence.
    %
    % REFERENCES:
    %  [1] D. T. Jamieson, and J. S. Tudhope, Desalination, 8, 393-401, 1970.
    %=========================================================================

    %% CHECK INPUT ARGUMENTS

    % CHECK THAT S&T HAVE SAME SHAPE
    if ~isequal(size(S),size(T))
        error('check_stp: S & T must have same dimensions');
    end

    % CONVERT TEMPERATURE INPUT TO °C
    switch lower(uT)
        case 'c'
        case 'k'
            T = T - 273.15;
        case 'f'
            T = 5/9*(T-32);
        case 'r'
            T = 5/9*(T-491.67);
        otherwise
            error('Not a recognized temperature unit. Please use ''C'', ''K'', ''F'', or ''R''');
    end

    % CONVERT SALINITY TO PPT
    switch lower(uS)
        case 'ppt'
        case 'ppm'
            S = S/1000;
        case 'w'
            S = S*1000;
        case '%'
            S = S*10;
        otherwise
            error('Not a recognized salinity unit. Please use ''ppt'', ''ppm'', ''w'', or ''%''');
    end

    % CHECK THAT S & T ARE WITHIN THE FUNCTION RANGE
    if ~isequal((T<0)+(T>180),zeros(size(T)))
        warning('Temperature is out of range for thermal conductivity function 0<T<180 C');
    end

    if ~isequal((S<0)+(S>160),zeros(size(S)))
        warning('Salinity is out of range for thermal conductivity function 0<S<160 g/kg');
    end

    %% BEGIN

    T = 1.00024*T;      %convert from T_90 to T_68
    S = S / 1.00472;    %convert from S to S_P
    k = 10.^(log10(240+0.0002*S)+0.434*(2.3-(343.5+0.037*S)./(T+273.15)).*(1-(T+273.15)./(647.3+0.03*S)).^(1/3)-3);

end
