classdef OrbitElement
    properties
        semimajor
        eccentricity
        inclination %Deg
        right_ascension %Deg
        argument_perigee %Deg
        true_anomaly %Deg
    end
    methods
        function obj = OrbitElement(a, e, i, Sigma, omega, theta)
            obj.semimajor = a;
            obj.eccentricity = e;
            obj.inclination = i;
            obj.right_ascension = Sigma;
            obj.argument_perigee = omega;
            obj.true_anomaly = theta;
        end

        function theta = trueanomalyAfter(obj, tp)
            arguments
                obj (1, 1) OrbitElement
                tp (1, 1) {mustBeReal, mustBeNonnegative}
            end

            global mu

            theta0 = deg2rad(obj.true_anomaly);
            a = obj.semimajor;
            e = obj.eccentricity;
            
            E0 = 2*atan(tan(theta0/2) * sqrt((1-e)/(1+e))); % Eccentric Anomaly
            M0 = (E0 - e*sin(E0)); % Mean Anomaly
            M1 = M0 + sqrt(mu/a^2) * tp;
            
            % function define for Newton-Raphson method
            f = @(E) M1 - (E - e*sin(E));
            df = @(E) -(1 - e*cos(E));
            E1 = NRmethod(E0, f, df);
            
            % true anoamly after tp
            theta = 2*atan(tan(E1/2) * sqrt((1+e)/(1-e))); 
        end

        function r = pos(obj, tp)
            arguments
                obj (1, 1) OrbitElement
                tp (1, 1) {mustBeReal, mustBeNonnegative}
            end

            a = obj.semimajor;
            e = obj.eccentricity;
            theta = obj.trueanomalyAfter(tp);
            r = a*(1-e^2)/(1+e*cos(theta)) * [cos(theta), sin(theta)];
        end

        function v = vel(obj, tp)
            arguments
                obj (1, 1) OrbitElement
                tp (1, 1) {mustBeReal, mustBeNonnegative}
            end

            global mu
            a = obj.semimajor;
            e = obj.eccentricity;
            theta = obj.trueanomalyAfter(tp);
            v = sqrt(mu/(a*(1-e^2))) * [-sin(theta), e+cos(theta)];
        end
    end
end