function result = lindistflow_down (casename, vmin, vmax, Ppcc)

%% ---------- Data Load ----------
    mpc = feval(casename);  
    define_constants; 
    baseMVA = mpc.baseMVA;

    bus = mpc.bus;
    branch = mpc.branch;
    branch = branch(branch(:, BR_STATUS) == 1, :); 
    nb = size(bus, 1);    
    nl = size(branch, 1);
    root = find(bus(:, BUS_TYPE)==3, 1); 

    Pd = bus(:, PD)/baseMVA;
    Qd = bus(:, QD)/baseMVA;
    r = branch(:, BR_R);
    x = branch(:, BR_X);
    Smax = branch(:, RATE_A)/baseMVA;
    Smax(~isfinite(Smax) | Smax <= 0) = inf;

    f = branch(:, F_BUS); t = branch(:, T_BUS);

    parentEdgeId = nan(nb, 1);
    for e = 1:nl
        parentEdgeId(t(e)) = e;  
    end

    PresMax = 2/100; 
    SresMax = 2/100;
    
    %% ---------- Variables ----------
    Pl = sdpvar(nl, 1); Ql = sdpvar(nl, 1); 
    p  = sdpvar(nb, 1); q  = sdpvar(nb, 1); 
    Pres = sdpvar(nb, 1); Qres = sdpvar(nb, 1);
    v2  = sdpvar(nb, 1); 

    %% ---------- Constraints ----------
    cons = [];
    
    % 1) nodal volt constraints
    cons = [cons; vmin^2 <= v2 <= vmax^2];
    %cons = [cons; v2(root) == bus(root, 8)^2]; 
    
    % 2) Res power outputs 
    N = 36;
    for n = 1:nb
        cons = [cons; 0<= Pres(n) <= PresMax];
        for k = 1:N
            cons = [cons; cos(k/N*pi)*Pres(n) + sin(k/N*pi)*Qres(n) <= SresMax];
            cons = [cons; -(cos(k/N*pi)*Pres(n) + sin(k/N*pi)*Qres(n)) <= SresMax];
        end
    end

    % 3) power balance constraints 
    for j = 1:nb
        if j ~= root
            if j ~= root + 1
                cons = [cons; p(j) == Pres(j)-Pd(j)];
                cons = [cons; q(j) == Qres(j)-Qd(j)];
            else
                cons = [cons; p(j) == 0];
                cons = [cons; q(j) == 0];
            end
        end
    end
    
    % 4) branch voltage drop constraints
    for e = 1:nl
        cons = [cons; v2(t(e)) == v2(f(e)) - 2*(r(e)*Pl(e) + x(e)*Ql(e))];
    end
    
    % 5) power flow constraints 
    
    % [Bus i] --Pij--> [Bus j] --Pjk1--> [Bus k1]
    %                    | \
    %                    |  \----Pjk2--> [Bus k2]
    %                    |
    %                    \-------Pjkn--> [Bus kn]
    %
    % LinDistFlow, ignore power loss
    
    for j = 1:nb
        if j == root % root node have no parent
            outEdges = find(f == j);
            cons = [cons; p(j)-sum(Pl(outEdges)) == 0];
            cons = [cons; q(j)-sum(Ql(outEdges)) == 0];
        else
            inEdge = parentEdgeId(j);
            outEdges = find(f == j);
            cons = [cons; p(j) + Pl(inEdge) - sum(Pl(outEdges)) == 0];
            cons = [cons; q(j) + Ql(inEdge) - sum(Ql(outEdges)) == 0];
        end
    end

    % 5) !! fix PCC active power for scanning method
    cons = [cons; Ppcc - 1e-5 <= p(root) <= Ppcc + 1e-5];

    %% solver settings
    ops = sdpsettings('solver', 'gurobi', 'verbose', 2);

    %% ---------- Objective ----------
    % minimize reactive power at PCC/root
    obj = q(root);
    sol = optimize(cons, obj, ops);

    if sol.problem ~= 0
        warning('Solver failed: %s', sol.info);
    end

    Times = sol.solvertime;
    
    %% ---------- Output ---------- 
    result.Pl    = value(Pl) * baseMVA;
    result.Ql    = value(Ql) * baseMVA;
    result.p     = value(p)  * baseMVA;
    result.q     = value(q)  * baseMVA;
    result.Pres  = value(Pres) * baseMVA;
    result.Qres  = value(Qres) * baseMVA;
    result.vn    = sqrt(value(v2));
    result.v2    = value(v2);
    result.obj   = value(obj) * baseMVA;
    
    result.Ppcc  = value(p(root)) * baseMVA;
    result.Qpcc  = value(q(root)) * baseMVA;
    result.Qmin  = value(q(root)) * baseMVA;
    result.status = sol.problem;
    result.info   = sol.info;
    result.time   = sol.solvertime;
    
end


