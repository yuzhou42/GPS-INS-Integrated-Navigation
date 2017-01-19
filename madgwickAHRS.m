  function  Quaternion = madgwickAHRS( Gyroscope, Accelerometer, Magnetometer, quaternion_l,Beta,SamplePeriod)
            q = quaternion_l; % short name local variable for readability
            % Normalise accelerometer measurement
            if(norm(Accelerometer) == 0), return; end	% handle NaN
            Accelerometer = Accelerometer / norm(Accelerometer);	% normalise magnitude

            % Normalise magnetometer measurement
            if(norm(Magnetometer) == 0), return; end	% handle NaN
            Magnetometer = Magnetometer / norm(Magnetometer);	% normalise magnitude

            % Reference direction of Earth's magnetic feild
            h = quaternProd(q, quaternProd([0 Magnetometer], quaternConj(q)));
            b = [0 norm([h(2) h(3)]) 0 h(4)];

            % Gradient decent algorithm corrective step
            F = [-2*(q(2)*q(4) - q(1)*q(3)) - Accelerometer(1)
                -2*(q(1)*q(2) + q(3)*q(4)) - Accelerometer(2)
                -2*(0.5 - q(2)^2 - q(3)^2) - Accelerometer(3)
                ((2*b(2)*(0.5 - q(3)^2 - q(4)^2) + 2*b(4)*(q(2)*q(4) - q(1)*q(3))) - Magnetometer(1))
                ((2*b(2)*(q(2)*q(3) - q(1)*q(4)) + 2*b(4)*(q(1)*q(2) + q(3)*q(4))) - Magnetometer(2))
                ((2*b(2)*(q(1)*q(3) + q(2)*q(4)) + 2*b(4)*(0.5 - q(2)^2 - q(3)^2)) - Magnetometer(3))];
            J = [2*q(3),                 	-2*q(4),                    2*q(1),                         -2*q(2)
                -2*q(2),                 	-2*q(1),                    	-2*q(4),                         -2*q(3)
                0,                         4*q(2),                    4*q(3),                         0
               -2*b(4)*q(3),               2*b(4)*q(4),               -4*b(2)*q(3)-2*b(4)*q(1),       -4*b(2)*q(4)+2*b(4)*q(2)
                -2*b(2)*q(4)+2*b(4)*q(2),	2*b(2)*q(3)+2*b(4)*q(1),	2*b(2)*q(2)+2*b(4)*q(4),       -2*b(2)*q(1)+2*b(4)*q(3)
                2*b(2)*q(3),                2*b(2)*q(4)-4*b(4)*q(2),	2*b(2)*q(1)-4*b(4)*q(3),        2*b(2)*q(2)];
            step = (J'*F);
            step = step / norm(step);	% normalise step magnitude

            % Compute rate of change of quaternion
            qDot = 0.5 * quaternProd(q, [0 Gyroscope]) - Beta * step';

            % Integrate to yield quaternion
            q = q + qDot * SamplePeriod;
            Quaternion = q / norm(q); % normalise quaternion          
        end