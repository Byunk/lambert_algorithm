classdef OrbitElement
    properties
        semimajor
        eccentricity
        inclination
        right_ascension
        argument_perigee
        true_anomaly
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
    end
end