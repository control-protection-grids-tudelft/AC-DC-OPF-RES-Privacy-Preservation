clc; clear; close all;

%% ---------- Settings ----------
casename = 'case32_res';
vmin = 0.95;
vmax = 1.05;

% Ppcc is in p.u. inside lindistflow_up/down
% Negative Ppcc means RES exports power from distribution network to PCC
Ppcc_set = 0:-0.05:-0.6;   % p.u., 

nP = length(Ppcc_set);

Ppcc_MW   = nan(nP,1);
Qpcc_max  = nan(nP,1);
Qpcc_min  = nan(nP,1);

%% ---------- P-sweep ----------
for k = 1:nP

    Ppcc = Ppcc_set(k);

    result_up   = lindistflow_up(casename, vmin, vmax, Ppcc);
    result_down = lindistflow_down(casename, vmin, vmax, Ppcc);

    if result_up.status == 0 && result_down.status == 0
        Ppcc_MW(k)  = result_up.Ppcc;
        Qpcc_max(k) = result_up.Qmax;
        Qpcc_min(k) = result_down.Qmin;
    else
        fprintf('Ppcc = %.4f p.u. infeasible or failed.\n', Ppcc);
    end
end

%% ---------- Remove infeasible points ----------
valid = ~isnan(Ppcc_MW) & ~isnan(Qpcc_max) & ~isnan(Qpcc_min);

Ppcc_MW  = Ppcc_MW(valid);
Qpcc_max = Qpcc_max(valid);
Qpcc_min = Qpcc_min(valid);

%% ---------- Plot export power as positive ----------
figure;
hold on; grid on; box on;

Pexport_MW = -Ppcc_MW;

P_env2 = [Pexport_MW; flipud(Pexport_MW)];
Q_env2 = [Qpcc_max; flipud(Qpcc_min)];

fill(P_env2, Q_env2, [0.85 0.85 0.85], ...
    'FaceAlpha', 0.4, ...
    'EdgeColor', 'none');

plot(Pexport_MW, Qpcc_max, 'o-', 'LineWidth', 1.5);
plot(Pexport_MW, Qpcc_min, 's-', 'LineWidth', 1.5);

xlabel('RES export power at PCC / MW');
ylabel('Q_{PCC} / MVAr');
legend('P-Q feasible region', 'Q_{PCC}^{max}', 'Q_{PCC}^{min}', ...
    'Location', 'best');

title('P-Q Feasible Region with Export Power Positive');


%% ---------- Convert PQ envelope to p.u. linear constraints ----------
mpc = feval(casename);
baseMVA = mpc.baseMVA;

% Convert MW / MVAr back to p.u.
Ppcc_pu   = Ppcc_MW   / baseMVA;
Qmax_pu   = Qpcc_max  / baseMVA;
Qmin_pu   = Qpcc_min  / baseMVA;

% Use Ppcc, not Pexport
points = [
    Ppcc_pu, Qmax_pu;
    flipud(Ppcc_pu), flipud(Qmin_pu)
];

% Remove duplicate points
points = unique(points, 'rows', 'stable');

% Convex hull
K = convhull(points(:,1), points(:,2));
vertices = points(K(1:end-1), :);

[A_pq, b_pq] = polygon_to_halfspace(vertices);

disp('A_pq = ');
disp(A_pq);

disp('b_pq = ');
disp(b_pq);

%% Example YALMIP usage
% Ppcc = sdpvar(1,1);
% Qpcc = sdpvar(1,1);
% cons = [cons; A_pq * [Ppcc; Qpcc] <= b_pq];


function [A_pq, b_pq] = polygon_to_halfspace(vertices)

    N = size(vertices, 1);
    
    % Ensure counter-clockwise order
    area = 0;
    for i = 1:N
        j = mod(i, N) + 1;
        area = area + vertices(i,1)*vertices(j,2) ...
                    - vertices(j,1)*vertices(i,2);
    end
    
    if area < 0
        vertices = flipud(vertices);
    end
    
    A_pq = zeros(N, 2);
    b_pq = zeros(N, 1);
    
    for i = 1:N
        j = mod(i, N) + 1;
    
        xi = vertices(i,1);
        yi = vertices(i,2);
        xj = vertices(j,1);
        yj = vertices(j,2);
    
        dx = xj - xi;
        dy = yj - yi;
    
        % For CCW polygon:
        % dy*x - dx*y <= dy*xi - dx*yi
        A_pq(i,:) = [dy, -dx];
        b_pq(i)   = dy*xi - dx*yi;
    end

end

