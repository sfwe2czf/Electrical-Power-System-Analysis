clc;


BUS = table2array(readtable("CASE_B_BUS DATA.csv"));
BRANCH = table2array(readtable("CASE_B_BRANCH DATA.csv"));
GEN = table2array(readtable("CASE_B_GEN DATA.csv"));

% Admittance Matrix
length_branch = length(BRANCH(:,1));
Number_bus  = length(BUS(:,1));
Number_gen = length(GEN(:,1));

Ybus = zeros(Number_bus, Number_bus);

for i = 1 : 1 : length_branch
    series_Z = BRANCH(i, 3) + 1i * BRANCH(i, 4);
    parallel_B = 1i * BRANCH(i, 5);

    from_bus = BRANCH(i, 1);
    to_bus = BRANCH(i, 2);

    Ybus(from_bus, from_bus) = Ybus(from_bus, from_bus) + 1/series_Z + parallel_B / 2;
    Ybus(from_bus, to_bus) = Ybus(from_bus, to_bus) - 1/series_Z;
    Ybus(to_bus, from_bus) = Ybus(to_bus, from_bus) - 1/series_Z;
    Ybus(to_bus, to_bus) = Ybus(to_bus, to_bus) + 1/series_Z + parallel_B/2;
end

for i = 1 : 1 : Number_bus
    Gs = BUS(i, 5)/100;
    Bs = BUS(i, 6)/100;

    Ybus(i, i) = Ybus(i, i) - Gs - 1i*Bs;
end

Gbus = real(Ybus);
Bbus = imag(Ybus);

% VARIABLE

Num_slack = 0;
Num_PV = 0;
Num_PQ = 0;

for i = 1 : 1 : Number_bus
    if (BUS(i, 2)== 3)
        Num_slack = Num_slack + 1;
    elseif (BUS(i, 2) == 2)
        Num_PV = Num_PV + 1;
    elseif (BUS(i, 2) == 1)
        Num_PQ = Num_PQ + 1;
    end
end

Theta = zeros(Num_PQ + Num_PV, 1);
Vmag = ones(Num_PQ, 1);
X = [Theta ; Vmag];

Theta_O = zeros(Number_bus, 1);
Vmag_O = ones(Number_bus, 1);

for i = 1 : 1 : Number_bus
    if (BUS(i, 2) == 3)
        Theta_O(i, 1) = BUS(i, 8);
    end
end

for i = 1 : 1 : Number_gen
    Vmag_O(GEN(i, 1), 1) = GEN(i, 4);
end

P_i_O = zeros (Number_bus, 1);
Q_i_O = zeros (Number_bus, 1);

for i = 1 : 1 : Number_gen
    P_i_O(GEN(i, 1), 1) = P_i_O(GEN(i, 1), 1) + GEN(i, 2)/100;
end

for i = 1 : 1 : Number_bus
    P_i_O(i, 1) = P_i_O(i, 1) - BUS(i, 3)/100;
    Q_i_O(i, 1) = Q_i_O(i, 1) - BUS(i, 4)/100;
end

Index_Bus = zeros(Num_PQ + Num_PV, 1);
Index_PV = 0;
Index_PQ = 0;

for i  = 1 : 1 : Number_bus
    if BUS(i, 2) == 3
    elseif BUS(i, 2) == 2
        Index_PV = Index_PV + 1;
        Index_Bus(Num_PQ + Index_PV, 1) = BUS(i, 1);
    elseif BUS(i, 2) == 1
        Index_PQ = Index_PQ + 1;
        Index_Bus(Index_PQ, 1) = BUS(i, 1);
    else
    end
end

P_i = zeros(Num_PQ + Num_PV, 1);
Q_i = zeros(Num_PQ, 1);

for i = 1 : 1 : Num_PQ + Num_PV
    P_i(i, 1) = P_i_O(Index_Bus(i, 1));
end

for i = 1 : 1 : Num_PQ
    Q_i(i, 1) = Q_i_O(Index_Bus(i, 1));
end

% Mismatch Vector

eps = 1;

while eps > 0.0001

    P_cal_O = zeros(Number_bus, 1);
    Q_cal_O = zeros(Number_bus, 1);

    for i = 1 : 1 : Num_PQ + Num_PV
        Theta(i, 1) = X(i, 1);
        Theta_O(Index_Bus(i, 1), 1) = X(i, 1);
    end

    for i = 1 : 1 : Num_PQ
        Vmag(i, 1) = X(Num_PQ+Num_PV+i, 1);
        Vmag_O(Index_Bus(i, 1), 1) = X(Num_PQ+Num_PV+i, 1);
    end

    for i = 1 : 1 : Number_bus
        for j = 1 : 1: Number_bus
            Theta_ij = Theta_O(i,1)-Theta_O(j, 1);
            P_cal_O(i, 1) = P_cal_O(i, 1) + Vmag_O(i, 1)*Vmag_O(j, 1)*(Gbus(i, j)*cos(Theta_ij) + Bbus(i, j)*sin(Theta_ij));
            Q_cal_O(i, 1) = Q_cal_O(i, 1) + Vmag_O(i, 1)*Vmag_O(j, 1)*(Gbus(i, j)*sin(Theta_ij) - Bbus(i, j)*cos(Theta_ij));
        end
    end
    
    P_cal = zeros(Num_PV + Num_PQ, 1);
    Q_cal = zeros(Num_PQ, 1);

    for i = 1 : 1 : Num_PQ + Num_PV
        P_cal(i, 1) = P_cal_O(Index_Bus(i, 1), 1);
        Theta(i, 1) = Theta_O(Index_Bus(i, 1), 1);
    end
    
    for i = 1 : 1 : Num_PQ
        Q_cal(i, 1) = Q_cal_O(Index_Bus(i, 1), 1);
        Vmag(i, 1) = Vmag_O(Index_Bus(i, 1), 1);
    end
    
    Delta_f = [P_cal - P_i ; Q_cal - Q_i];

    Jacobian_P_Theta = zeros(Num_PQ + Num_PV, Num_PQ + Num_PV);
    Jacobian_P_V = zeros(Num_PQ + Num_PV, Num_PQ);
    Jacobian_Q_Theta = zeros(Num_PQ, Num_PQ + Num_PV);
    Jacobian_Q_V = zeros(Num_PQ, Num_PQ);

    % for Jacobian_P_Theta
    for i = 1 : 1 : Num_PV + Num_PQ
        for j = 1 : 1 : Num_PQ + Num_PV
            Theta_ij = Theta(i, 1) - Theta(j, 1);
            if i == j
                Jacobian_P_Theta(i, j) = -Q_cal_O(Index_Bus(i, 1), 1) - Bbus(Index_Bus(i, 1), Index_Bus(j, 1)) * Vmag_O(Index_Bus(i, 1), 1)^2;
            else
                Jacobian_P_Theta(i, j) = Vmag_O(Index_Bus(i, 1), 1)*Vmag_O(Index_Bus(j, 1), 1)*(Gbus(Index_Bus(i, 1), Index_Bus(j, 1))*sin(Theta_ij) - Bbus(Index_Bus(i, 1), Index_Bus(j, 1))*cos(Theta_ij));
            end
        end
    end
    
    % for Jacobian_P_V
    for i = 1 : 1 : Num_PV + Num_PQ
        for j = 1 : 1 : Num_PQ
            Theta_ij = Theta(i, 1) - Theta(j, 1);
            if i == j
                Jacobian_P_V(i, j) = P_cal_O(Index_Bus(i, 1), 1)/Vmag_O(Index_Bus(i, 1), 1) + Gbus(Index_Bus(i, 1), Index_Bus(j, 1)) * Vmag_O(Index_Bus(i, 1), 1);
            else
                Jacobian_P_V(i, j) = Vmag_O(Index_Bus(i, 1), 1)*(Gbus(Index_Bus(i, 1), Index_Bus(j, 1))*cos(Theta_ij) + Bbus(Index_Bus(i, 1), Index_Bus(j, 1))*sin(Theta_ij));
            end
        end
    end

    % for Jacobian_Q_Theta
    for i = 1 : 1 : Num_PQ
        for j = 1 : 1 : Num_PQ + Num_PV
            Theta_ij = Theta(i, 1) - Theta(j, 1);
            if i == j
                Jacobian_Q_Theta(i, j) = P_cal_O(Index_Bus(i, 1), 1) - Gbus(Index_Bus(i, 1), Index_Bus(j, 1)) * Vmag_O(Index_Bus(i, 1), 1)^2;
            else
                Jacobian_Q_Theta(i, j) = -Vmag_O(Index_Bus(i, 1), 1)*Vmag_O(Index_Bus(j, 1), 1)*(Gbus(Index_Bus(i, 1), Index_Bus(j, 1))*cos(Theta_ij) + Bbus(Index_Bus(i, 1), Index_Bus(j, 1))*sin(Theta_ij));
            end
        end
    end

    % for Jacobian_Q_V
    for i = 1 : 1 : Num_PQ
        for j = 1 : 1 : Num_PQ
            Theta_ij = Theta(i, 1) - Theta(j, 1);
            if i == j
                Jacobian_Q_V(i, j) = Q_cal_O(Index_Bus(i, 1), 1)/Vmag_O(Index_Bus(i, 1), 1) - Bbus(Index_Bus(i, 1), Index_Bus(j, 1)) * Vmag_O(Index_Bus(i, 1), 1);
            else
                Jacobian_Q_V(i, j) = Vmag_O(Index_Bus(i, 1), 1)*(Gbus(Index_Bus(i, 1), Index_Bus(j, 1))*sin(Theta_ij) - Bbus(Index_Bus(i, 1), Index_Bus(j, 1))*cos(Theta_ij));
            end
        end
    end

    Jacobian_Full = [Jacobian_P_Theta, Jacobian_P_V ; Jacobian_Q_Theta, Jacobian_Q_V];

    Delta_X = -inv(Jacobian_Full) * Delta_f;

    X = X + Delta_X;
    eps = max(abs(Delta_X));

end

Result = [BUS(:, 1), Vmag_O, Theta_O*180/pi];
disp(Result);
