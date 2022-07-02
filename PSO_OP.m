function out = PSO_OP(problem, parameters, data, size_init)
% minimizing cost function for optimal parameters in FLC
% Problem contains PSO algorithm parameters.
% parameters contains fuzzy logic controller(FLC) parameters.
% data contains main project details.
% updated in ver. 18 based on parEPSO ver. 17

disp('PSO_OP running...')

%% Initiallization
rule_set_zahedi = data.rule_set_zahedi;
OR.optimal_rule_set = data.rule_set_zahedi;

%% Definition
load = problem.load;         % load vevtor in kW
nVar = problem.param_nVar;   % search space dimension of FLC parameters
Var_size = [1 nVar];         % matrix size of decision variables

% lower bound of decision variables
Var_min = problem.parameters_min;
% upper bound of decision variables
Var_max = problem.parameters_max;

%% Parameters
max_iteration = parameters.max_iteration; % max number of iterations
nPop = parameters.nPop;         % population size

w_min = parameters.w_min;       % min inertia coefficient
w_max = parameters.w_max;       % max inertia coefficient
Kw1 = parameters.Kw1;            % shape-shift coefficient for w for x
Kw2 = parameters.Kw2;            % shape-shift coefficient for w for parameters
Kw3 = parameters.Kw3;            % shape-shift coefficient for w for rules

c1_min = parameters.c1_min;     % min personal acceleration coefficient
c1_max = parameters.c1_max;     % max personal acceleration coefficient
Kc1 = parameters.Kc1;           % shape-shift coefficient for c1

c2_min = parameters.c2_min;     % min social acceleration coefficient
c2_max = parameters.c2_max;     % max social acceleration coefficient
Kc2 = parameters.Kc2;           % shape-shift coefficient for c2

Max_velocity = 0.2*(Var_max - Var_min);  % max velocity
Min_velocity = - Max_velocity;           % min velocity

% flag for showing best solution in each iteration
show_iteration = parameters.show_iteration;

%% Initialization
% flag for stucking in local extremun
stuck = 0;
% itteration number that PSO breaks / for callback
break_it = 0;

% a particle template
empty_particle.position = [];
empty_particle.velocity = [];
empty_particle.cost = [];
empty_particle.best.position = [];
empty_particle.best.cost = [];

% create population array
particle = repmat(empty_particle, nPop, 1);
% PB is the buffer for particle.best during parallel computation
PB = particle.best;

% initial global best
global_best.cost = inf;

% detemine size for parfor loop
global_best.position = zeros(Var_size);

% best_cost vector
particles_cost = inf*ones(1,nPop);

% parpool('local')
% initial poulation members
parfor i=1:nPop
    particle(i).position = [sort(rand(1,7)), sort(rand(1,10)),...
        sort(rand(1,7)), sort(rand(1,10)), sort(rand(1,10)), sort(rand(1,7)),...
        sort(rand(1,10)), sort(rand(1,10)), sort(rand(1,7))];
    
    % initial velocity
    particle(i).velocity = zeros(Var_size);
    
    % evaluation
    x = [size_init, particle(i).position, rule_set_zahedi];
    particle(i).cost = CostFunction(x,data);
        
    % initial update personal best
    particle(i).best.position = particle(i).position;
    particle(i).best.cost = particle(i).cost;
end

% updating global best cost and position in parallel
tmp = [particle.best];
[PB, GB_index] = sort([tmp.cost],'ascend');
global_best = tmp(GB_index(1));

% array to hold the best cost function value at each iteration
best_costs = zeros(max_iteration, 1);

%% Main Loop of PSO
for it=1:max_iteration
    fprintf('PSO_OP it #%d running\n',it);
    w = 0.5*(w_max+w_min) + Kw3*(atan( pi + it*(-2*pi/max_iteration))...
        ) * (w_max-w_min);
    c1 = 0.5*(c1_max+c1_min) + Kc1*(atan( pi + it*(-2*pi/max_iteration))...
        ) * (c1_max-c1_min);
    c2 = 0.5*(c2_max+c2_min) + Kc2*(atan( -pi + it*(2*pi/max_iteration))...
        ) * (c2_max-c2_min);
    
    % applying w limits
    w = max(w_min, w);
    w = min(w_max, w);
    
    % applying c1 limits
    c1 = max(c1_min, c1);
    c1 = min(c1_max, c1);
    
    % applying c2 limits
    c2 = max(c2_min, c2);
    c2 = min(c2_max, c2);
        
    tic
    parfor i=1:nPop
        % update velocity
        particle(i).velocity = w*particle(i).velocity ...
            + c1*rand(Var_size).*(particle(i).best.position - particle(i).position) ...
            + c2*rand(Var_size).*(global_best.position - particle(i).position);
        
        % apply velocity limits
        particle(i).velocity = max(particle(i).velocity, Min_velocity);
        particle(i).velocity = min(particle(i).velocity, Max_velocity);
        
        % update position
        particle(i).position = particle(i).position + particle(i).velocity;
        
        % apply lower bound and upper ound limits
        particle(i).position = max(particle(i).position, Var_min);
        particle(i).position = min(particle(i).position, Var_max);
        
        % sorting particle(i).position for FLC parameters
        particle(i).position = SortPosition(particle(i).position);
                
        % evaluation
        x = [size_init, particle(i).position, rule_set_zahedi];
        [particle(i).cost, detail] = CostFunction(x,data);
        particles_cost(i) = particle(i).cost;
        
        % solution constraint check
        particle(i).flag = Constraints(x,data,detail);
        
        % update personal best
        if (particle(i).cost < particle(i).best.cost) && particle(i).flag
            particle(i).best.position = particle(i).position;
            particle(i).best.cost = particle(i).cost;
        end
    end
    toc
    
    % update personal best
    parfor i=1:nPop
        if particle(i).cost<particle(i).best.cost
            particle(i).best.position = particle(i).position;
            particle(i).best.cost = particle(i).cost;
        end
    end

    % update global best
        % updating global best position in parallel
        % if min(particles_cost)<GB
        %     global_best.cost = min(particles_cost);
        %     GB = global_best.cost;
        % end
    tmp = [particle.best];
    [PB, GB_index] = sort([tmp.cost],'ascend');
        % global_best = tmp(GB_index(1));

    if (PB(1)<global_best.cost) && particle(GB_index(1)).flag
        global_best = tmp(GB_index(1));
        % calculate optimal rule_set for the new global optimum
        OR.stuck=1;
        try_it = 1;
        % optimal rule_set
        while (OR.stuck==1) && (try_it<parameters.try_it)
            disp('GB improved');
            OR = parPSO_OR(problem, parameters, data, [size_init, global_best.position]);
            try_it = try_it + 1;
        end
    end
    
    % store the best cost of the value
    best_costs(it) = global_best.cost;
    
%     % display iteration information
%     if show_iteration
%         disp(['Iteration ' num2str(it) ...
%             ': Best cost = ' num2str(best_costs(it)) ...
%             ' w = ' num2str(w) ' c1 = ' num2str(c1) ...
%             ' c2 = ' num2str(c2)]);
%     end
    
        if it>parameters.stuck_step
            % checking the improvement
            if best_costs(it) == best_costs(it-parameters.stuck_step)
                break_it = it;
                % checking fot stucking in local extremum
                if best_costs(it) == best_costs(1)
                    stuck = 1;
                end
                break
            end
        end
end

% checking if PSO is stuck in a local extermum
if best_costs(1)==best_costs(end)
    stuck=1;
end

out.population = particle;
out.best_solution = global_best;
out.optimal_parameters = global_best.position;
out.optimal_rule_set = OR.optimal_rule_set;
out.best_costs = best_costs;
out.stuck = stuck;
out.break_it = break_it;

%% Logging best_cost for paper (with total cost(TC))
if data.excel==1
    filename = 'details_FLC_log.xlsx';
    sheet = 1;
    xlRange = 'A1';
    
    % loading the excel file to retrieve previous results
    last_file = xlsread(filename);
    
    % changing the first value of the last log(written by CostFunction.m)
    full = size(last_file,2);
    blanked = full - length(out.best_costs) -length(out.optimal_parameters) -length(out.optimal_rule_set) -6;
    new_file = [last_file;
                1111*ones(1,full);
                out.stuck,data.w,parameters.max_iteration,...
                parameters.nPop,parameters.stuck_step,parameters.try_it,...
                out.best_costs',zeros(1,blanked),out.optimal_parameters,...
                out.optimal_rule_set;
                1111*ones(1,full)];
    xlswrite(filename,new_file,sheet,xlRange);
end

disp('PSO_OP ended.')
end

