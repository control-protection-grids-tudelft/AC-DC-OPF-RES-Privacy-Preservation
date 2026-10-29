clc; clear; close all;

%% S1.1 Load AC network parameters & variables 

% -------------------------------------------------------------------------
% 1.1.1 AC grid parameters
% -------------------------------------------------------------------------
mpc_ac = case9;                       
bus_ac     = mpc_ac.bus;
gen_ac     = mpc_ac.gen;
branch_ac  = mpc_ac.branch;
baseMVA_ac = mpc_ac.baseMVA;
gencost_ac = mpc_ac.gencost;

[Ybus_ac, ~, ~] = makeYbus(baseMVA_ac, bus_ac, branch_ac);
Ybus_ac = full(Ybus_ac);

G_ac = real(Ybus_ac);
B_ac = imag(Ybus_ac);

nbus_ac    = size(bus_ac, 1);
ngen_ac    = size(gen_ac, 1);
nbranch_ac = size(branch_ac, 1);

pd_ac = bus_ac(:, 3) / baseMVA_ac;
qd_ac = bus_ac(:, 4) / baseMVA_ac;

refBus_ac = find(bus_ac(:, 2) == 3);    
genBus_ac = gen_ac(:, 1);  
a2vBus_ac = [4; 6; 8];

% -------------------------------------------------------------------------
% 1.1.2 AC grid variables
% -------------------------------------------------------------------------
pn_ac       = sdpvar(nbus_ac, 1);
qn_ac       = sdpvar(nbus_ac, 1);
pgen_ac     = sdpvar(nbus_ac, 1);
qgen_ac     = sdpvar(nbus_ac, 1);
pa2v_ac     = sdpvar(nbus_ac, 1);
qa2v_ac     = sdpvar(nbus_ac, 1);
pl_ac       = sdpvar(nbus_ac, nbus_ac, 'full');
ql_ac       = sdpvar(nbus_ac, nbus_ac, 'full');
cc_ac       = sdpvar(nbus_ac, nbus_ac, 'full');
ss_ac       = sdpvar(nbus_ac, nbus_ac, 'full');


%% S1.2 Add SOCP-relaxed AC network constraints

cons_ac = [];

% -------------------------------------------------------------------------
% 1.2.1 Branch active/reactive power flow equations
%   Smax: apparent power limits derived from MATPOWER branch ratings (rateA)
%   Polygonal approximation is used to enforce |S_ij| <= Smax(i,j),
%   i.e., p_ij^2 + q_ij^2 <= Smax(i,j)^2
% -------------------------------------------------------------------------
Smax_ac = zeros(nbus_ac, nbus_ac);

for k = 1:nbranch_ac
    i = branch_ac(k, 1);
    j = branch_ac(k, 2);

    rateA = branch_ac(k, 6);

    if rateA == 0
        Smax_val = 1e3;   
    else
        Smax_val = rateA / baseMVA_ac;
    end

    Smax_ac(i, j) = Smax_val;
    Smax_ac(j, i) = Smax_val;   
end

N = 36;
for k = 1:nbranch_ac
    i = branch_ac(k, 1);
    j = branch_ac(k, 2);

    cons_ac = [cons_ac, pl_ac(i, j) == -G_ac(i, j)*(cc_ac(i, i) ...
        -cc_ac(i, j))-B_ac(i, j)*ss_ac(i, j)];
    cons_ac = [cons_ac, ql_ac(i, j) == B_ac(i, j)*(cc_ac(i, i) ...
        -cc_ac(i, j))-G_ac(i, j)*ss_ac(i, j)];
    % ij and ji
    cons_ac = [cons_ac, pl_ac(j, i) == -G_ac(j, i)*(cc_ac(j, j) ...
        -cc_ac(j, i))-B_ac(j, i)*ss_ac(j, i)];
    cons_ac = [cons_ac, ql_ac(j, i) == B_ac(j, i)*(cc_ac(j, j) ...
        -cc_ac(j, i))-G_ac(j, i)*ss_ac(j, i)];
    
    for k = 1: N
        cons_ac = [cons_ac, cos(k/N*pi)*pl_ac(i, j) ...
            + sin(k/N*pi)*ql_ac(i, j) <= Smax_ac(i, j)];
        cons_ac = [cons_ac; -(cos(k/N*pi)*pl_ac(i, j) ...
            + sin(k/N*pi)*ql_ac(i, j)) <= Smax_ac(i, j)];
        % ij and ji
        cons_ac = [cons_ac, cos(k/N*pi)*pl_ac(j, i) ...
            + sin(k/N*pi)*ql_ac(j, i) <= Smax_ac(j, i)];
        cons_ac = [cons_ac; -(cos(k/N*pi)*pl_ac(j, i) ...
            + sin(k/N*pi)*ql_ac(j, i)) <= Smax_ac(j, i)];
    end
    
end

% -------------------------------------------------------------------------
% 1.2.2 Nodal active/reactive power flow equations
% -------------------------------------------------------------------------
for i = 1:nbus_ac
    cons_ac = [cons_ac, pn_ac(i) == ss_ac(i, i)*B_ac(i, i) ...
        + sum(cc_ac(i, :).*G_ac(i, :) - ss_ac(i, :).*B_ac(i, :))];
    cons_ac = [cons_ac, qn_ac(i) == ss_ac(i, i)*G_ac(i, i) ...
        - sum(cc_ac(i, :).*B_ac(i, :) + ss_ac(i, :).*G_ac(i, :))];
end

% -------------------------------------------------------------------------
% 1.2.3 Fix injections of all non-slack buses
%     IMPORTANT:
%     Slack bus injections are NOT fixed, so that the slack bus can balance
%     power mismatch/losses, which is consistent with standard PF logic.
% -------------------------------------------------------------------------
for i = 1:nbus_ac
    % if ~ismember(i, refBus_ac)
        p_expr = -pd_ac(i);
        q_expr = -qd_ac(i);

        if ismember(i, genBus_ac)
            p_expr = p_expr + pgen_ac(i);
            q_expr = q_expr + qgen_ac(i);
        end

        if ismember(i, a2vBus_ac)
            p_expr = p_expr - pa2v_ac(i);
            q_expr = q_expr - qa2v_ac(i);
        end

        cons_ac = [cons_ac, pn_ac(i) == p_expr];
        cons_ac = [cons_ac, qn_ac(i) == q_expr];

    % end
end

% -------------------------------------------------------------------------
% 1.2.4 SOCP relaxation constraints
%     Symmetry constraints
% -------------------------------------------------------------------------
for k = 1:nbranch_ac
    i = branch_ac(k, 1);
    j = branch_ac(k, 2);
    cons_ac = [cons_ac, cc_ac(i, j)^2 + ss_ac(i, j)^2 ...
        <= cc_ac(i, i)*cc_ac(j, j)];
end

for i = 1:nbus_ac
    cons_ac = [cons_ac, ss_ac(i, i) == 0];  
    cons_ac=[cons_ac, cc_ac(i, i)>=0];
    for j = i + 1:nbus_ac
        cons_ac = [cons_ac, cc_ac(i, j) == cc_ac(j, i)];
        cons_ac = [cons_ac, ss_ac(i, j) + ss_ac(j, i) == 0];
    end
end


% -------------------------------------------------------------------------
% 1.2.5 Generator active/reactive power limits
% -------------------------------------------------------------------------
pgenMin_ac = zeros(nbus_ac, 1);
pgenMax_ac = zeros(nbus_ac, 1);
qgenMin_ac = zeros(nbus_ac, 1);
qgenMax_ac = zeros(nbus_ac, 1);

for i = 1:nbus_ac
    gen_idx = find(gen_ac(:,1) == i, 1, 'first');
    if ~isempty(gen_idx)
        qgenMax_ac(i) = gen_ac(gen_idx, 4)  / baseMVA_ac;   
        qgenMin_ac(i) = gen_ac(gen_idx, 5)  / baseMVA_ac;   
        pgenMax_ac(i) = gen_ac(gen_idx, 9)  / baseMVA_ac;   
        pgenMin_ac(i) = gen_ac(gen_idx, 10) / baseMVA_ac;   
    end
end

for i = 1:nbus_ac
    if ismember(i, genBus_ac)
        cons_ac = [cons_ac, pgenMin_ac(i) <= pgen_ac(i) <= pgenMax_ac(i)];
        cons_ac = [cons_ac, qgenMin_ac(i) <= qgen_ac(i) <= qgenMax_ac(i)];
    else
        cons_ac = [cons_ac, pgen_ac(i) == 0];
        cons_ac = [cons_ac, qgen_ac(i) == 0];
    end
end

% -------------------------------------------------------------------------
% 1.2.6 Voltage magnitude bounds
% -------------------------------------------------------------------------
for i = 1:nbus_ac
    Vmax_ac = bus_ac(i, 12);
    Vmin_ac = bus_ac(i, 13);
    cons_ac = [cons_ac, Vmin_ac^2 <= cc_ac(i, i) <= Vmax_ac^2];
end

% -------------------------------------------------------------------------
% 1.2.7 Voltage references
% -------------------------------------------------------------------------
r = refBus_ac.'
gen_idx = find(gen_ac(:,1) == r, 1, 'first');
Vref_ac = gen_ac(gen_idx,6);
cons_ac = [cons_ac, cc_ac(r, r) == Vref_ac^2];

% -------------------------------------------------------------------------
% 1.2.8 Obj: Minizing generation costs
% -------------------------------------------------------------------------

c2_ac = zeros(nbus_ac, 1);
c1_ac = zeros(nbus_ac, 1);
c0_ac = zeros(nbus_ac, 1);

for i = 1:nbus_ac
    gen_idx = find(gen_ac(:,1) == i, 1, 'first');
    if ~isempty(gen_idx)
        ncost = gencost_ac(gen_idx, 4);

        if ncost == 3
            c2_ac(i) = gencost_ac(gen_idx, end-2);
            c1_ac(i) = gencost_ac(gen_idx, end-1);
            c0_ac(i) = gencost_ac(gen_idx, end);
        elseif ncost == 2
            c2_ac(i) = 0;
            c1_ac(i) = gencost_ac(gen_idx, end-1);
            c0_ac(i) = gencost_ac(gen_idx, end);
        end
    end
end

obj_ac = 0;

for i = 1:nbus_ac
    if ismember(i, genBus_ac)
        Pg_MW = pgen_ac(i) * 1;
        obj_ac = obj_ac + c2_ac(i) * Pg_MW^2 + c1_ac(i) * Pg_MW ...
            + c0_ac(i)/baseMVA_ac;
    end
end

% -------------------------------------------------------------------------
% 1.2.9 Solve
% -------------------------------------------------------------------------
% opt = sdpsettings('solver','gurobi', ...
%                   'verbose',1, ...
%                   'gurobi.QCPDual',1);
% 
% sol_ac = optimize(cons_ac, obj_ac, opt);
% 
% if sol_ac.problem ~= 0
%     disp('Optimization failed.');
%     disp(sol_ac.info);
%     yalmiperror(sol_ac.problem);
%     return;
% end
% 
% disp('AC Optimization succeeded.');


%% S2.1 Load DC network parameters & variables 

% -------------------------------------------------------------------------
% 2.1.1 DC grid parameters
% -------------------------------------------------------------------------
mpc_dc      = case5_dc;                       
results_dc  = runpf(mpc_dc);
bus_dc      = results_dc.bus;
gen_dc      = results_dc.gen;
branch_dc   = results_dc.branch;
baseMW_dc   = results_dc.baseMVA;

[Ybus_dc, ~, ~] = makeYbus(baseMW_dc, bus_dc, branch_dc);
Ybus_dc = full(Ybus_dc);

nbus_dc         = size(bus_dc, 1);
ngen_dc         = size(gen_dc, 1);
nbranch_dc      = size(branch_dc, 1);

% Construct equivalent "impedance"
zij_dc = 1 ./ abs(Ybus_dc);
zij_dc = zij_dc - diag(diag(zij_dc));
zij_dc(isinf(zij_dc)) = 1e5;

bus_dc(:, 3) = bus_dc(:, 3) / baseMW_dc; 
pd_dc = bus_dc(:, 3);

fbus_dc = branch_dc(:, 1);
tbus_dc = branch_dc(:, 2);

refBus_dc = find(bus_dc(:, 2) == 3);     

% -------------------------------------------------------------------------
% 2.1.2 DC grid parameters
% -------------------------------------------------------------------------
u_dc        = sdpvar(nbus_dc, 1);              
lij_dc      = sdpvar(nbus_dc, nbus_dc, 'full');   
pn_dc       = sdpvar(nbus_dc, 1);              
pij_dc      = sdpvar(nbus_dc, nbus_dc, 'full');   

%% S2.2 Add SOCP-relaxed AC network constraints
cons_dc=[]; 
M = 1e5;

% -------------------------------------------------------------------------
% 2.2.1 Nodal power balance
% -------------------------------------------------------------------------
for i = 1:nbus_dc
    expr = 0;
    for k = 1:nbranch_dc
        if fbus_dc(k) == i
            j = tbus_dc(k);
            expr = expr + pij_dc(i, j);
        elseif tbus_dc(k) == i
            j = fbus_dc(k);
            expr = expr + pij_dc(i, j);
        end
    end
    cons_dc = [cons_dc; pn_dc(i) == expr];
end

% -------------------------------------------------------------------------
% 2.2.2 Branch equations
% -------------------------------------------------------------------------
for k = 1:nbranch_dc
    i = fbus_dc(k);
    j = tbus_dc(k);
    rateA = branch_dc(k, 6);

    if rateA == 0
        Pmax_val = 1e3;   
    else
        Pmax_val = rateA / baseMW_dc;
    end
    Pmax_val = 100 / baseMW_dc;

    % loss
    cons_dc = [cons_dc;
        pij_dc(i, j) + pij_dc(j, i) == zij_dc(i, j) * lij_dc(i, j)];

    % volt drop
    cons_dc = [cons_dc;
       u_dc(i) - u_dc(j) == zij_dc(i, j) * (pij_dc(i, j) - pij_dc(j, i))];

    % SOCP
    cons_dc = [cons_dc;
        pij_dc(i, j)^2 <= lij_dc(i, j) * u_dc(i);
        pij_dc(j, i)^2 <= lij_dc(j, i) * u_dc(j)];

    cons_dc = [cons_dc;
        lij_dc(i, j) >= 0
        lij_dc(j, i) >= 0];

    cons_dc = [cons_dc; 
        u_dc(i) >= 0];

    cons_dc = [cons_dc;
        lij_dc(j, i) == lij_dc(i, j)];

    % Bounds
    cons_dc = [cons_dc;
       pij_dc(i, j) <= Pmax_val
       pij_dc(j, i) <= Pmax_val];

end

% ---------------------------------------------------------------------
% 2.2.3 Power Injection
% -------------------------------------------------------------------------
% for i = 1:nbus_dc
%     if ~ismember(i, refBus_dc)
%         cons_dc = [cons_dc; pn_dc(i) == -pd_dc(i)];
%     end
% end

% -------------------------------------------------------------------------
% 2.2.4 Slack bus voltage
% -------------------------------------------------------------------------
% r = refBus_dc.'
% gen_idx = find(gen_dc(:,1) == r, 1, 'first');
% Vref_dc = gen_dc(gen_idx,6);
% cons_dc = [cons_dc, u_dc(r) == Vref_dc^2];

% -------------------------------------------------------------------------
% 2.2.5 Voltage magnitude bounds
% -------------------------------------------------------------------------
for i = 1:nbus_dc
    Vmax_dc = bus_dc(i,12);
    Vmin_dc = bus_dc(i,13);
    cons_dc = [cons_dc, Vmin_dc^2 <= u_dc(i) <= Vmax_dc^2];
end


%% S3.1 RES1 parameters & variables 
cons_res1 = [];
% -------------------------------------------------------------------------
% 3.1.1 RES parameters 
% -------------------------------------------------------------------------
mpc_res1     = case32_res            
define_constants;
baseMVA_res1 = mpc_res1.baseMVA;

bus_res1     = mpc_res1.bus;
branch_res1  = mpc_res1.branch;
branch_res1  = branch_res1(branch_res1(:, BR_STATUS) == 1, :);   

nb_res1      = size(bus_res1, 1);
nl_res1      = size(branch_res1, 1);
root_res1    = find(bus_res1(:, BUS_TYPE) == 3, 1);

r_res1       = branch_res1(:, BR_R);
x_res1       = branch_res1(:, BR_X);

Pd_res1 = bus_res1(:, PD)/baseMVA_res1;
Qd_res1 = bus_res1(:, QD)/baseMVA_res1;
Smax_res1    = branch_res1(:, RATE_A) / baseMVA_res1;
Smax_res1(~isfinite(Smax_res1) | Smax_res1 <= 0) = inf;

f_res1       = branch_res1(:, F_BUS);
t_res1       = branch_res1(:, T_BUS);

parentEdgeId_res1 = nan(nb_res1, 1);
for e = 1:nl_res1
    parentEdgeId_res1(t_res1(e)) = e;  
end

PresMax_res1 = 2/100; 
SresMax_res1 = 2/100;

vmin_res1 = 0.95; 
vmax_res1 = 1.05;

%% ---------- Variables ----------
Pl_res1 = sdpvar(nl_res1, 1);     Ql_res1 = sdpvar(nl_res1, 1); 
p_res1  = sdpvar(nb_res1, 1);     q_res1  = sdpvar(nb_res1, 1); 
Pres_res1 = sdpvar(nb_res1, 1);   Qres_res1 = sdpvar(nb_res1, 1);
v2_res1  = sdpvar(nb_res1, 1); 

pr2v_res1    = sdpvar(1, 1);
qr2v_res1    = sdpvar(1, 1);

%% ---------- Constraints ----------
cons_res1 = [];
    
% 1) nodal volt constraints
cons_res1 = [cons_res1; vmin_res1^2 <= v2_res1 <= vmax_res1^2];
%cons_res1 = [cons_res1; v2_res1(root_res1) == bus_res1(root_res1, 8)^2]; 
    
% 2) Res power outputs 
N = 36;
for n = 1:nb_res1
    cons_res1 = [cons_res1; 0<= Pres_res1(n) <= PresMax_res1];
    for k = 1:N
        cons_res1 = [cons_res1; cos(k/N*pi)*Pres_res1(n) ...
            + sin(k/N*pi)*Qres_res1(n) <= SresMax_res1];
        cons_res1 = [cons_res1; -(cos(k/N*pi)*Pres_res1(n) ...
            + sin(k/N*pi)*Qres_res1(n)) <= SresMax_res1];
    end
end

% 3) power balance constraints 
for j = 1:nb_res1
    if j ~= root_res1
        if j ~= root_res1 + 1
            cons_res1 = [cons_res1; p_res1(j) == Pres_res1(j)-Pd_res1(j)];
            cons_res1 = [cons_res1; q_res1(j) == Qres_res1(j)-Qd_res1(j)];
        else
            cons_res1 = [cons_res1; p_res1(j) == 0];
            cons_res1 = [cons_res1; q_res1(j) == 0]; 
        end
    end
 end
    
% 4) branch voltage drop constraints
for e = 1:nl_res1
    cons_res1 = [cons_res1; v2_res1(t_res1(e)) == v2_res1(f_res1(e)) ...
    - 2*(r_res1(e)*Pl_res1(e) + x_res1(e)*Ql_res1(e))];
end
    
% 5) power flow constraints 
    
% [Bus i] --Pij--> [Bus j] --Pjk1--> [Bus k1]
%                    | \
%                    |  \----Pjk2--> [Bus k2]
%                    |
%                    \-------Pjkn--> [Bus kn]
%
% LinDistFlow, ignore power loss
    
for j = 1:nb_res1
    if j == root_res1 % root node have no parent
        outEdges = find(f_res1 == j);
        cons_res1 = [cons_res1; p_res1(j)-sum(Pl_res1(outEdges)) == 0];
        cons_res1 = [cons_res1; q_res1(j)-sum(Ql_res1(outEdges)) == 0];
    else
        inEdge = parentEdgeId_res1(j);
        outEdges = find(f_res1 == j);
        cons_res1 = [cons_res1; p_res1(j) + Pl_res1(inEdge) ...
                - sum(Pl_res1(outEdges)) == 0];
        cons_res1 = [cons_res1; q_res1(j) + Ql_res1(inEdge) ...
                - sum(Ql_res1(outEdges)) == 0];
    end
end

% 6) power coupling constraints 
cons_res1 = [cons_res1; p_res1(root_res1) == -pr2v_res1];
cons_res1 = [cons_res1; q_res1(root_res1) == -qr2v_res1];


%% S3.2 RES2 parameters & variables 
cons_res2 = [];
% -------------------------------------------------------------------------
% 3.2.1 RES parameters 
% -------------------------------------------------------------------------
mpc_res2     = case32_res            
define_constants;
baseMVA_res2 = mpc_res2.baseMVA;

bus_res2     = mpc_res2.bus;
branch_res2  = mpc_res2.branch;
branch_res2  = branch_res2(branch_res2(:, BR_STATUS) == 1, :);   

nb_res2      = size(bus_res2, 1);
nl_res2      = size(branch_res2, 1);
root_res2    = find(bus_res2(:, BUS_TYPE) == 3, 1);

r_res2       = branch_res2(:, BR_R);
x_res2       = branch_res2(:, BR_X);

Pd_res2 = bus_res2(:, PD)/baseMVA_res2;
Qd_res2 = bus_res2(:, QD)/baseMVA_res2;
Smax_res2    = branch_res2(:, RATE_A) / baseMVA_res2;
Smax_res2(~isfinite(Smax_res2) | Smax_res2 <= 0) = inf;

f_res2       = branch_res2(:, F_BUS);
t_res2       = branch_res2(:, T_BUS);

parentEdgeId_res2 = nan(nb_res2, 1);
for e = 1:nl_res2
    parentEdgeId_res2(t_res2(e)) = e;  
end

PresMax_res2 = 2/100; 
SresMax_res2 = 2/100;

vmin_res2 = 0.95; 
vmax_res2 = 1.05;

%% ---------- Variables ----------
Pl_res2 = sdpvar(nl_res2, 1);     Ql_res2 = sdpvar(nl_res2, 1); 
p_res2  = sdpvar(nb_res2, 1);     q_res2  = sdpvar(nb_res2, 1); 
Pres_res2 = sdpvar(nb_res2, 1);   Qres_res2 = sdpvar(nb_res2, 1);
v2_res2  = sdpvar(nb_res2, 1); 

pr2v_res2    = sdpvar(1, 1);
qr2v_res2    = sdpvar(1, 1);

%% ---------- Constraints ----------
cons_res2 = [];
    
% 1) nodal volt constraints
cons_res2 = [cons_res2; vmin_res2^2 <= v2_res2 <= vmax_res2^2];
%cons_res2 = [cons_res2; v2_res2(root_res2) == bus_res2(root_res2, 8)^2]; 
    
% 2) Res power outputs 
N = 36;
for n = 1:nb_res2
    cons_res2 = [cons_res2; 0<= Pres_res2(n) <= PresMax_res2];
    for k = 1:N
        cons_res2 = [cons_res2; cos(k/N*pi)*Pres_res2(n) ...
            + sin(k/N*pi)*Qres_res2(n) <= SresMax_res2];
        cons_res2 = [cons_res2; -(cos(k/N*pi)*Pres_res2(n) ...
            + sin(k/N*pi)*Qres_res2(n)) <= SresMax_res2];
    end
end

% 3) power balance constraints 
for j = 1:nb_res2
    if j ~= root_res2
        if j ~= root_res2 + 1
            cons_res2 = [cons_res2; p_res2(j) == Pres_res2(j)-Pd_res2(j)];
            cons_res2 = [cons_res2; q_res2(j) == Qres_res2(j)-Qd_res2(j)];
        else
            cons_res2 = [cons_res2; p_res2(j) == 0];
            cons_res2 = [cons_res2; q_res2(j) == 0]; 
        end
    end
 end
    
% 4) branch voltage drop constraints
for e = 1:nl_res2
    cons_res2 = [cons_res2; v2_res2(t_res2(e)) == v2_res2(f_res2(e)) ...
    - 2*(r_res2(e)*Pl_res2(e) + x_res2(e)*Ql_res2(e))];
end
    
% 5) power flow constraints 
    
% [Bus i] --Pij--> [Bus j] --Pjk1--> [Bus k1]
%                    | \
%                    |  \----Pjk2--> [Bus k2]
%                    |
%                    \-------Pjkn--> [Bus kn]
%
% LinDistFlow, ignore power loss
    
for j = 1:nb_res2
    if j == root_res2 % root node have no parent
        outEdges = find(f_res2 == j);
        cons_res2 = [cons_res2; p_res2(j)-sum(Pl_res2(outEdges)) == 0];
        cons_res2 = [cons_res2; q_res2(j)-sum(Ql_res2(outEdges)) == 0];
    else
        inEdge = parentEdgeId_res2(j);
        outEdges = find(f_res2 == j);
        cons_res2 = [cons_res2; p_res2(j) + Pl_res2(inEdge) ...
                - sum(Pl_res2(outEdges)) == 0];
        cons_res2 = [cons_res2; q_res2(j) + Ql_res2(inEdge) ...
                - sum(Ql_res2(outEdges)) == 0];
    end
end

% 6) power coupling constraints 
cons_res2 = [cons_res2; p_res2(root_res2) == -pr2v_res2];
cons_res2 = [cons_res2; q_res2(root_res2) == -qr2v_res2];


%% S4.1 Load VSC network parameters & variables 

% -------------------------------------------------------------------------
% 4.1.1 VSC parameters
% -------------------------------------------------------------------------
mpc_vsc1 = case2_vsc;                       
bus_vsc1     = mpc_vsc1.bus;
branch_vsc1  = mpc_vsc1.branch;
baseMVA_vsc1 = mpc_vsc1.baseMVA;

[Ybus_vsc1, ~, ~] = makeYbus(baseMVA_vsc1, bus_vsc1, branch_vsc1);
Ybus_vsc1 = full(Ybus_vsc1);

G_vsc1 = real(Ybus_vsc1);
B_vsc1 = imag(Ybus_vsc1);

nbus_vsc1    = size(bus_vsc1, 1);
nbranch_vsc1 = size(branch_vsc1, 1);

pd_vsc1 = bus_vsc1(:, 3) / baseMVA_vsc1;
qd_vsc1 = bus_vsc1(:, 4) / baseMVA_vsc1;

refBus_vsc1 = find(bus_vsc1(:, 2) == 3);    


mpc_vsc2 = case2_vsc;                       
bus_vsc2     = mpc_vsc2.bus;
branch_vsc2  = mpc_vsc2.branch;
baseMVA_vsc2 = mpc_vsc2.baseMVA;

[Ybus_vsc2, ~, ~] = makeYbus(baseMVA_vsc2, bus_vsc2, branch_vsc2);
Ybus_vsc2 = full(Ybus_vsc2);

G_vsc2 = real(Ybus_vsc2);
B_vsc2 = imag(Ybus_vsc2);

nbus_vsc2    = size(bus_vsc2, 1);
nbranch_vsc2 = size(branch_vsc2, 1);

pd_vsc2 = bus_vsc2(:, 3) / baseMVA_vsc2;
qd_vsc2 = bus_vsc2(:, 4) / baseMVA_vsc2;

refBus_vsc2 = find(bus_vsc2(:, 2) == 3);   


mpc_vsc3 = case2_vsc;                       
bus_vsc3     = mpc_vsc3.bus;
branch_vsc3  = mpc_vsc3.branch;
baseMVA_vsc3 = mpc_vsc3.baseMVA;

[Ybus_vsc3, ~, ~] = makeYbus(baseMVA_vsc3, bus_vsc3, branch_vsc3);
Ybus_vsc3 = full(Ybus_vsc3);

G_vsc3 = real(Ybus_vsc3);
B_vsc3 = imag(Ybus_vsc3);

nbus_vsc3    = size(bus_vsc3, 1);
nbranch_vsc3 = size(branch_vsc3, 1);

pd_vsc3 = bus_vsc2(:, 3) / baseMVA_vsc3;
qd_vsc3 = bus_vsc2(:, 4) / baseMVA_vsc3;

refBus_vsc3 = find(bus_vsc3(:, 2) == 3);


mpc_vsc4 = case2_vsc;                       
bus_vsc4     = mpc_vsc4.bus;
branch_vsc4  = mpc_vsc4.branch;
baseMVA_vsc4 = mpc_vsc4.baseMVA;

[Ybus_vsc4, ~, ~] = makeYbus(baseMVA_vsc4, bus_vsc4, branch_vsc4);
Ybus_vsc4 = full(Ybus_vsc4);

G_vsc4 = real(Ybus_vsc4);
B_vsc4 = imag(Ybus_vsc4);

nbus_vsc4    = size(bus_vsc4, 1);
nbranch_vsc4 = size(branch_vsc4, 1);

pd_vsc4 = bus_vsc2(:, 3) / baseMVA_vsc4;
qd_vsc4 = bus_vsc2(:, 4) / baseMVA_vsc4;

refBus_vsc4 = find(bus_vsc4(:, 2) == 3);


mpc_vsc5 = case2_vsc;                       
bus_vsc5     = mpc_vsc5.bus;
branch_vsc5  = mpc_vsc5.branch;
baseMVA_vsc5 = mpc_vsc5.baseMVA;

[Ybus_vsc5, ~, ~] = makeYbus(baseMVA_vsc5, bus_vsc5, branch_vsc5);
Ybus_vsc5 = full(Ybus_vsc5);

G_vsc5 = real(Ybus_vsc5);
B_vsc5 = imag(Ybus_vsc5);

nbus_vsc5    = size(bus_vsc5, 1);
nbranch_vsc5 = size(branch_vsc5, 1);

pd_vsc5 = bus_vsc2(:, 3) / baseMVA_vsc5;
qd_vsc5 = bus_vsc2(:, 4) / baseMVA_vsc5;

refBus_vsc5 = find(bus_vsc5(:, 2) == 3);

aloss_vsc=0.0011; bloss_vsc=0.003; closs_vsc=0.0046;
Imax_vsc = 1;
sMax_vsc1 = 100 /baseMVA_vsc1 ;
sMax_vsc2 = 100 /baseMVA_vsc2 ;
sMax_vsc3 = 100 /baseMVA_vsc3 ;
sMax_vsc4 = 100 /baseMVA_vsc4 ;
sMax_vsc5 = 100 /baseMVA_vsc5 ;

% -------------------------------------------------------------------------
% 4.1.2 VSC variables
% -------------------------------------------------------------------------
pn_vsc1       = sdpvar(nbus_vsc1, 1);
qn_vsc1       = sdpvar(nbus_vsc1, 1);
pgen_vsc1     = sdpvar(nbus_vsc1, 1);
qgen_vsc1     = sdpvar(nbus_vsc1, 1);
pa2v_vsc1     = sdpvar(nbus_vsc1, 1);
qa2v_vsc1     = sdpvar(nbus_vsc1, 1);
pl_vsc1       = sdpvar(nbus_vsc1, nbus_vsc1, 'full');
ql_vsc1       = sdpvar(nbus_vsc1, nbus_vsc1, 'full');
cc_vsc1       = sdpvar(nbus_vsc1, nbus_vsc1, 'full');
ss_vsc1       = sdpvar(nbus_vsc1, nbus_vsc1, 'full');
pv2d_vsc1     = sdpvar(1, 1);
pv2a_vsc1     = sdpvar(1, 1);
qv2a_vsc1     = sdpvar(1, 1);
I_vsc1        = sdpvar(1, 1);
l_vsc1        = sdpvar(1, 1);
ploss_vsc1    = sdpvar(1, 1);


pn_vsc2       = sdpvar(nbus_vsc2, 1);
qn_vsc2       = sdpvar(nbus_vsc2, 1);
pgen_vsc2     = sdpvar(nbus_vsc2, 1);
qgen_vsc2     = sdpvar(nbus_vsc2, 1);
pa2v_vsc2     = sdpvar(nbus_vsc2, 1);
qa2v_vsc2     = sdpvar(nbus_vsc2, 1);
pl_vsc2       = sdpvar(nbus_vsc2, nbus_vsc2, 'full');
ql_vsc2       = sdpvar(nbus_vsc2, nbus_vsc2, 'full');
cc_vsc2       = sdpvar(nbus_vsc2, nbus_vsc2, 'full');
ss_vsc2       = sdpvar(nbus_vsc2, nbus_vsc2, 'full');
pv2d_vsc2     = sdpvar(1, 1);
pv2a_vsc2     = sdpvar(1, 1);
qv2a_vsc2     = sdpvar(1, 1);
I_vsc2        = sdpvar(1, 1);
l_vsc2        = sdpvar(1, 1);
ploss_vsc2    = sdpvar(1, 1);


pn_vsc3       = sdpvar(nbus_vsc3, 1);
qn_vsc3       = sdpvar(nbus_vsc3, 1);
pgen_vsc3     = sdpvar(nbus_vsc3, 1);
qgen_vsc3     = sdpvar(nbus_vsc3, 1);
pa2v_vsc3     = sdpvar(nbus_vsc3, 1);
qa2v_vsc3     = sdpvar(nbus_vsc3, 1);
pl_vsc3       = sdpvar(nbus_vsc3, nbus_vsc3, 'full');
ql_vsc3       = sdpvar(nbus_vsc3, nbus_vsc3, 'full');
cc_vsc3       = sdpvar(nbus_vsc3, nbus_vsc3, 'full');
ss_vsc3       = sdpvar(nbus_vsc3, nbus_vsc3, 'full');
pv2d_vsc3     = sdpvar(1, 1);
pv2a_vsc3     = sdpvar(1, 1);
qv2a_vsc3     = sdpvar(1, 1);
I_vsc3        = sdpvar(1, 1);
l_vsc3        = sdpvar(1, 1);
ploss_vsc3    = sdpvar(1, 1);


pn_vsc4       = sdpvar(nbus_vsc4, 1);
qn_vsc4       = sdpvar(nbus_vsc4, 1);
pgen_vsc4     = sdpvar(nbus_vsc4, 1);
qgen_vsc4     = sdpvar(nbus_vsc4, 1);
pa2v_vsc4     = sdpvar(nbus_vsc4, 1);
qa2v_vsc4     = sdpvar(nbus_vsc4, 1);
pl_vsc4       = sdpvar(nbus_vsc4, nbus_vsc4, 'full');
ql_vsc4       = sdpvar(nbus_vsc4, nbus_vsc4, 'full');
cc_vsc4       = sdpvar(nbus_vsc4, nbus_vsc4, 'full');
ss_vsc4       = sdpvar(nbus_vsc4, nbus_vsc4, 'full');
pv2d_vsc4     = sdpvar(1, 1);
pv2a_vsc4     = sdpvar(1, 1);
qv2a_vsc4     = sdpvar(1, 1);
I_vsc4        = sdpvar(1, 1);
l_vsc4        = sdpvar(1, 1);
ploss_vsc4    = sdpvar(1, 1);


pn_vsc5       = sdpvar(nbus_vsc5, 1);
qn_vsc5       = sdpvar(nbus_vsc5, 1);
pgen_vsc5     = sdpvar(nbus_vsc5, 1);
qgen_vsc5     = sdpvar(nbus_vsc5, 1);
pa2v_vsc5     = sdpvar(nbus_vsc5, 1);
qa2v_vsc5     = sdpvar(nbus_vsc5, 1);
pl_vsc5       = sdpvar(nbus_vsc5, nbus_vsc5, 'full');
ql_vsc5       = sdpvar(nbus_vsc5, nbus_vsc5, 'full');
cc_vsc5       = sdpvar(nbus_vsc5, nbus_vsc5, 'full');
ss_vsc5       = sdpvar(nbus_vsc5, nbus_vsc5, 'full');
pv2d_vsc5     = sdpvar(1, 1);
pv2a_vsc5     = sdpvar(1, 1);
qv2a_vsc5     = sdpvar(1, 1);
I_vsc5        = sdpvar(1, 1);
l_vsc5        = sdpvar(1, 1);
ploss_vsc5    = sdpvar(1, 1);


%% S4.2 VSC constraints

cons_vsc1 = [];
cons_vsc2 = [];
cons_vsc3 = [];
cons_vsc4 = [];
cons_vsc5 = [];

% -------------------------------------------------------------------------
% 4.2.1 Nodal active/reactive power flow equations
% -------------------------------------------------------------------------
for i = 1:nbus_vsc1
    cons_vsc1 = [cons_vsc1; pn_vsc1(i) == ss_vsc1(i, i)*B_vsc1(i, i) ...
        + sum(cc_vsc1(i, :).*G_vsc1(i, :) - ss_vsc1(i, :).*B_vsc1(i, :))];
    cons_vsc1 = [cons_vsc1; qn_vsc1(i) == ss_vsc1(i, i)*G_vsc1(i, i) ...
        - sum(cc_vsc1(i, :).*B_vsc1(i, :) + ss_vsc1(i, :).*G_vsc1(i, :))];
    cons_vsc1 = [cons_vsc1; .9^2 <= cc_vsc1(i, i) <= 1.1^2];
end

for i = 1:nbus_vsc2
    cons_vsc2 = [cons_vsc2; pn_vsc2(i) == ss_vsc2(i, i)*B_vsc2(i, i) ...
        + sum(cc_vsc2(i, :).*G_vsc2(i, :) - ss_vsc2(i, :).*B_vsc2(i, :))];
    cons_vsc2 = [cons_vsc2; qn_vsc2(i) == ss_vsc2(i, i)*G_vsc2(i, i) ...
        - sum(cc_vsc2(i, :).*B_vsc2(i, :) + ss_vsc2(i, :).*G_vsc2(i, :))];
    cons_vsc2 = [cons_vsc2; .9^2 <= cc_vsc2(i, i) <= 1.1^2];
end

for i = 1:nbus_vsc3
    cons_vsc3 = [cons_vsc3; pn_vsc3(i) == ss_vsc3(i, i)*B_vsc3(i, i) ...
        + sum(cc_vsc3(i, :).*G_vsc3(i, :) - ss_vsc3(i, :).*B_vsc3(i, :))];
    cons_vsc3 = [cons_vsc3; qn_vsc3(i) == ss_vsc3(i, i)*G_vsc3(i, i) ...
        - sum(cc_vsc3(i, :).*B_vsc3(i, :) + ss_vsc3(i, :).*G_vsc3(i, :))];
    cons_vsc3 = [cons_vsc3; .9^2 <= cc_vsc3(i, i) <= 1.1^2];
end

for i = 1:nbus_vsc4
    cons_vsc4 = [cons_vsc4; pn_vsc4(i) == ss_vsc4(i, i)*B_vsc4(i, i) ...
        + sum(cc_vsc4(i, :).*G_vsc4(i, :) - ss_vsc4(i, :).*B_vsc4(i, :))];
    cons_vsc4 = [cons_vsc4; qn_vsc4(i) == ss_vsc4(i, i)*G_vsc4(i, i) ...
        - sum(cc_vsc4(i, :).*B_vsc4(i, :) + ss_vsc4(i, :).*G_vsc4(i, :))];
    cons_vsc4 = [cons_vsc4; .9^2 <= cc_vsc4(i, i) <= 1.1^2];
end

for i = 1:nbus_vsc5
    cons_vsc5 = [cons_vsc5; pn_vsc5(i) == ss_vsc5(i, i)*B_vsc5(i, i) ...
        + sum(cc_vsc5(i, :).*G_vsc5(i, :) - ss_vsc5(i, :).*B_vsc5(i, :))];
    cons_vsc5 = [cons_vsc5; qn_vsc5(i) == ss_vsc5(i, i)*G_vsc5(i, i) ...
        - sum(cc_vsc5(i, :).*B_vsc5(i, :) + ss_vsc5(i, :).*G_vsc5(i, :))];
    cons_vsc5 = [cons_vsc5; .9^2 <= cc_vsc5(i, i) <= 1.1^2];
end


% -------------------------------------------------------------------------
% 4.2.2 PCC side power balance
% -------------------------------------------------------------------------

cons_vsc1 = [cons_vsc1; pn_vsc1(1) == -pv2a_vsc1];
cons_vsc1 = [cons_vsc1; qn_vsc1(1) == -qv2a_vsc1];    

cons_vsc2 = [cons_vsc2; pn_vsc2(1) == -pv2a_vsc2];
cons_vsc2 = [cons_vsc2; qn_vsc2(1) == -qv2a_vsc2];    

cons_vsc3 = [cons_vsc3; pn_vsc3(1) == -pv2a_vsc3];
cons_vsc3 = [cons_vsc3; qn_vsc3(1) == -qv2a_vsc3];    

cons_vsc4 = [cons_vsc4; pn_vsc4(1) == -pv2a_vsc4];
cons_vsc4 = [cons_vsc4; qn_vsc4(1) == -qv2a_vsc4];    

cons_vsc5 = [cons_vsc5; pn_vsc5(1) == -pv2a_vsc5];
cons_vsc5 = [cons_vsc5; qn_vsc5(1) == -qv2a_vsc5];    

% -------------------------------------------------------------------------
% 4.2.3 Capacity cosntraint
% -------------------------------------------------------------------------

N = 36;

for k = 1: N
        cons_vsc1 = [cons_vsc1; cos(k/N*pi)*pn_vsc1(2) ...
            + sin(k/N*pi)*qn_vsc1(2) <= sMax_vsc1];
        cons_vsc1 = [cons_vsc1; -(cos(k/N*pi)*pn_vsc1(2) ...
            + sin(k/N*pi)*qn_vsc1(2)) <= sMax_vsc1];

        cons_vsc2 = [cons_vsc2; cos(k/N*pi)*pn_vsc2(2) ...
            + sin(k/N*pi)*qn_vsc2(2) <= sMax_vsc2];
        cons_vsc2 = [cons_vsc2; -(cos(k/N*pi)*pn_vsc2(2) ...
            + sin(k/N*pi)*qn_vsc2(2)) <= sMax_vsc2];


        cons_vsc3 = [cons_vsc3; cos(k/N*pi)*pn_vsc3(2) ...
            + sin(k/N*pi)*qn_vsc3(2) <= sMax_vsc3];
        cons_vsc3 = [cons_vsc3; -(cos(k/N*pi)*pn_vsc3(2) ...
            + sin(k/N*pi)*qn_vsc3(2)) <= sMax_vsc3];

        cons_vsc4 = [cons_vsc4; cos(k/N*pi)*pn_vsc4(2) ...
            + sin(k/N*pi)*qn_vsc4(2) <= sMax_vsc4];
        cons_vsc4 = [cons_vsc4; -(cos(k/N*pi)*pn_vsc4(2) ...
            + sin(k/N*pi)*qn_vsc4(2)) <= sMax_vsc4];

        cons_vsc5 = [cons_vsc5; cos(k/N*pi)*pn_vsc5(2) ...
            + sin(k/N*pi)*qn_vsc5(2) <= sMax_vsc5];
        cons_vsc5 = [cons_vsc5; -(cos(k/N*pi)*pn_vsc5(2) ...
            + sin(k/N*pi)*qn_vsc5(2)) <= sMax_vsc5];
end


% -------------------------------------------------------------------------
% 4.2.4 VSC power loss
% -------------------------------------------------------------------------
cons_vsc1 = [cons_vsc1; pn_vsc1(2) + ploss_vsc1 + pv2d_vsc1 == 0;
    ploss_vsc1 == aloss_vsc + bloss_vsc*I_vsc1 +closs_vsc*l_vsc1;
    0 <= I_vsc1 <= Imax_vsc;
    0 <= l_vsc1 <= Imax_vsc^2
    ];

cons_vsc2 = [cons_vsc2; pn_vsc2(2) + ploss_vsc2 + pv2d_vsc2 == 0;
    ploss_vsc2 == aloss_vsc + bloss_vsc*I_vsc2 +closs_vsc*l_vsc2;
    0 <= I_vsc2 <= Imax_vsc;
    0 <= l_vsc2 <= Imax_vsc^2
 ];

cons_vsc3 = [cons_vsc3; pn_vsc3(2) + ploss_vsc3 + pv2d_vsc3 == 0;
    ploss_vsc3 == aloss_vsc + bloss_vsc*I_vsc3 +closs_vsc*l_vsc3;
    0 <= I_vsc3 <= Imax_vsc;
    0 <= l_vsc3 <= Imax_vsc^2
 ];
  
cons_vsc4 = [cons_vsc4; pn_vsc4(2) + ploss_vsc4 + pv2d_vsc4 == 0;
    ploss_vsc4 == aloss_vsc + bloss_vsc*I_vsc4 +closs_vsc*l_vsc4;
    0 <= I_vsc4 <= Imax_vsc;
    0 <= l_vsc4 <= Imax_vsc^2
 ];

cons_vsc5 = [cons_vsc5; pn_vsc5(2) + ploss_vsc5 + pv2d_vsc5 == 0;
    ploss_vsc5 == aloss_vsc + bloss_vsc*I_vsc5 +closs_vsc*l_vsc5;
    0 <= I_vsc5 <= Imax_vsc;
    0 <= l_vsc5 <= Imax_vsc^2
 ];


% -------------------------------------------------------------------------
% 4.2.5 SOCP relaxation
% -------------------------------------------------------------------------
cons_vsc1 = [cons_vsc1; I_vsc1^2 <= l_vsc1;
   pn_vsc1(2)^2 + qn_vsc1(2)^2 <= cc_vsc1(2, 2) * l_vsc1];

for k = 1:nbranch_vsc1
    i = branch_vsc1(k, 1);
    j = branch_vsc1(k, 2);
    cons_vsc1 = [cons_vsc1, cc_vsc1(i, j)^2 + ss_vsc1(i, j)^2 ...
        <= cc_vsc1(i, i)*cc_vsc1(j, j)];
end

for i = 1:nbus_vsc1
    cons_vsc1 = [cons_vsc1, ss_vsc1(i, i) == 0];  
    cons_vsc1=[cons_vsc1, cc_vsc1(i, i)>=0];
    for j = i + 1:nbus_vsc1
        cons_vsc1 = [cons_vsc1, cc_vsc1(i, j) == cc_vsc1(j, i)];
        cons_vsc1 = [cons_vsc1, ss_vsc1(i, j) + ss_vsc1(j, i) == 0];
    end
end

cons_vsc2 = [cons_vsc2; I_vsc2^2 <= l_vsc2;
   pn_vsc2(2)^2 + qn_vsc2(2)^2 <= cc_vsc2(2, 2) * l_vsc2];

for k = 1:nbranch_vsc2
    i = branch_vsc2(k, 1);
    j = branch_vsc2(k, 2);
    cons_vsc2 = [cons_vsc2, cc_vsc2(i, j)^2 + ss_vsc2(i, j)^2 ...
        <= cc_vsc2(i, i)*cc_vsc2(j, j)];
end

for i = 1:nbus_vsc2
    cons_vsc2 = [cons_vsc2, ss_vsc2(i, i) == 0];  
    cons_vsc2=[cons_vsc2, cc_vsc2(i, i)>=0];
    for j = i + 1:nbus_vsc2
        cons_vsc2 = [cons_vsc2, cc_vsc2(i, j) == cc_vsc2(j, i)];
        cons_vsc2 = [cons_vsc2, ss_vsc2(i, j) + ss_vsc2(j, i) == 0];
    end
end

cons_vsc3 = [cons_vsc3; I_vsc3^2 <= l_vsc3;
   pn_vsc3(2)^2 + qn_vsc3(2)^2 <= cc_vsc3(2, 2) * l_vsc3];

for k = 1:nbranch_vsc3
    i = branch_vsc3(k, 1);
    j = branch_vsc3(k, 2);
    cons_vsc3 = [cons_vsc3, cc_vsc3(i, j)^2 + ss_vsc3(i, j)^2 ...
        <= cc_vsc3(i, i)*cc_vsc3(j, j)];
end

for i = 1:nbus_vsc3
    cons_vsc3 = [cons_vsc3, ss_vsc3(i, i) == 0];  
    cons_vsc3=[cons_vsc3, cc_vsc3(i, i)>=0];
    for j = i + 1:nbus_vsc3
        cons_vsc3 = [cons_vsc3, cc_vsc3(i, j) == cc_vsc3(j, i)];
        cons_vsc3 = [cons_vsc3, ss_vsc3(i, j) + ss_vsc3(j, i) == 0];
    end
end

cons_vsc4 = [cons_vsc4; I_vsc4^2 <= l_vsc4;
   pn_vsc4(2)^2 + qn_vsc4(2)^2 <= cc_vsc4(2, 2) * l_vsc4];

for k = 1:nbranch_vsc4
    i = branch_vsc4(k, 1);
    j = branch_vsc4(k, 2);
    cons_vsc4 = [cons_vsc4, cc_vsc4(i, j)^2 + ss_vsc4(i, j)^2 ...
        <= cc_vsc4(i, i)*cc_vsc4(j, j)];
end

for i = 1:nbus_vsc4
    cons_vsc4 = [cons_vsc4, ss_vsc4(i, i) == 0];  
    cons_vsc4 =[cons_vsc4, cc_vsc4(i, i)>=0];
    for j = i + 1:nbus_vsc4
        cons_vsc4 = [cons_vsc4, cc_vsc4(i, j) == cc_vsc4(j, i)];
        cons_vsc4 = [cons_vsc4, ss_vsc4(i, j) + ss_vsc4(j, i) == 0];
    end
end

cons_vsc5 = [cons_vsc5; I_vsc5^2 <= l_vsc5;
   pn_vsc5(2)^2 + qn_vsc5(2)^2 <= cc_vsc5(2, 2) * l_vsc5];

for k = 1:nbranch_vsc5
    i = branch_vsc5(k, 1);
    j = branch_vsc5(k, 2);
    cons_vsc5 = [cons_vsc5, cc_vsc5(i, j)^2 + ss_vsc5(i, j)^2 ...
        <= cc_vsc5(i, i)*cc_vsc5(j, j)];
end

for i = 1:nbus_vsc5
    cons_vsc5 = [cons_vsc5, ss_vsc5(i, i) == 0];  
    cons_vsc5=[cons_vsc5, cc_vsc5(i, i)>=0];
    for j = i + 1:nbus_vsc5
        cons_vsc5 = [cons_vsc5, cc_vsc5(i, j) == cc_vsc5(j, i)];
        cons_vsc5 = [cons_vsc5, ss_vsc5(i, j) + ss_vsc5(j, i) == 0];
    end
end

% -------------------------------------------------------------------------
% 4.2.6 Obj
% -------------------------------------------------------------------------
obj_vsc = sum(pn_vsc1.^2) + sum(pn_vsc2.^2) + sum(pn_vsc3.^2) ...
    + sum(pn_vsc4.^2) + sum(pn_vsc5.^2)


%% S5 Coupling constraints
cons_cp = [];

cons_cp = [cons_cp; pa2v_ac(4) + pv2a_vsc1 == 0;
     qa2v_ac(4) + qv2a_vsc1 == 0;
     cc_ac(4, 4) == cc_vsc1(1, 1)];

cons_cp = [cons_cp; pa2v_ac(6) + pv2a_vsc2 == 0;
     qa2v_ac(6) + qv2a_vsc2 == 0;
     cc_ac(6, 6) == cc_vsc2(1, 1)];

cons_cp = [cons_cp; pa2v_ac(8) + pv2a_vsc3 == 0;
     qa2v_ac(8) + qv2a_vsc3 == 0;
     cc_ac(8, 8) == cc_vsc3(1, 1)];

cons_cp = [cons_cp; pr2v_res1 + pv2a_vsc4 == 0;
     qr2v_res1 + qv2a_vsc4 == 0;
     v2_res1(1) == cc_vsc4(1, 1)];
 
cons_cp = [cons_cp; pr2v_res2 + pv2a_vsc5 == 0;
     qr2v_res2 + qv2a_vsc5 == 0;
     v2_res2(1) == cc_vsc5(1, 1)];

cons_cp = [cons_cp; pn_dc(1) - pv2d_vsc1 == 0;
    pn_dc(2) - pv2d_vsc2 == 0;
    pn_dc(3) - pv2d_vsc3 == 0;
    pn_dc(4) - pv2d_vsc4 == 0;
    pn_dc(5) - pv2d_vsc5 == 0];



%% S6 Merge
cons_vsc     = [cons_vsc1; cons_vsc2; cons_vsc3; cons_vsc4; cons_vsc5];
cons_res     = [cons_res1; cons_res2];
cons         = [cons_ac; cons_dc; cons_vsc; cons_res; cons_cp];
obj          = obj_ac ;

opt = sdpsettings('solver', 'gurobi', 'verbose',1);
sol = optimize(cons, obj, opt);

if sol.problem ~= 0
    disp('Optimization failed.');
    disp(sol.info);
    yalmiperror(sol.problem);
    return;
end

value(obj)
solvertime = sol.solvertime

Pg_MW = value(pgen_ac(genBus_ac)) * baseMVA_ac;
generation_cost = sum(c2_ac(genBus_ac) .* Pg_MW.^2 ...
    + c1_ac(genBus_ac) .* Pg_MW + c0_ac(genBus_ac));
fprintf('\nTotal generation cost ($/h): %.6f\n', generation_cost);

disp(' ')
disp('================ AC Generator Outputs (MW) ================')

Pg = value(pgen_ac) * baseMVA_ac;

for i = 1:nbus_ac
    if ismember(i, genBus_ac)
        fprintf('Generator at Bus %2d : %8.3f MW\n', ...
            i, Pg(i));
    end
end


disp(' ')
disp('================ DC Port Powers (MW) =====================')

Pdc = value(pn_dc) * baseMW_dc;

for i = 1:nbus_dc
    fprintf('DC Bus %2d : %8.3f MW\n', ...
        i, Pdc(i));
end



