function out = parEPSO_OP(problem, parameters, data, size_init)
% out = parEPSO_OP(problem, parameters, data, size_init)
% minimizing cost function for optimal parameters in FLC
% Problem contains PSO algorithm parameters.
% parameters contains fuzzy logic controller(FLC) parameters.
% data contains main project details.

disp('parEPSO_OP running...')

%% Initiallization
rule_set_initial = data.rule_set_initial;
OR.optimal_rule_set = data.rule_set_initial;
parameters_initial = data.parameters_initial;

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

CR = parameters.CR;             % crossover rate (probability)
F_min = parameters.F_min;        % minimum random factor
F_max = parameters.F_max;       % maximum random factor

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

% initial mutual vector
M = empty_particle;

% initial trial vector
T = empty_particle;

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
parfor k=1:nPop
    particle(k).position = [sort(rand(1,7)), sort(rand(1,10)),...
        sort(rand(1,7)), sort(rand(1,10)), sort(rand(1,10)), sort(rand(1,7)),...
        sort(rand(1,10)), sort(rand(1,10)), sort(rand(1,7))];
    
    % initial velocity
    particle(k).velocity = zeros(Var_size);
    
    % evaluation
    x = [size_init, particle(k).position, rule_set_initial];
    particle(k).cost = CostFunction(x,data);
    
    % initial update personal best
    particle(k).best.position = particle(k).position;
    particle(k).best.cost = particle(k).cost;
end

% using the un-optimized parameters as initial data (particle(1))
if parameters.use_initial_guess
    particle(1).position = parameters_initial;
    
    % evaluation
    x = [size_init, parameters_initial, rule_set_initial];
    particle(1).cost = CostFunction(x,data);
        
    % initial update personal best
    particle(1).best.position = parameters_initial;
    particle(1).best.cost = particle(1).cost;
end

% updating global best cost and position in parallel
tmp = [particle.best];
[PB, GB_index] = sort([tmp.cost],'ascend');
global_best = tmp(GB_index(1));

% F, r, and b_rand are used in making mutual and trial vector
F = F_min + (F_max - F_min) * rand();  % scaling factor

r = randperm(nPop, 3);                 % random numbers(r1, r2, and r3)

%% making mutual vector and trial vector
% mutual vector M
rnd = rand();
M.position = Var_min + rnd()*(Var_max - Var_min);
% apply lower bound and upper ound limits
M.position = max(M.position, 0.99*Var_min);
M.position = min(M.position, 0.99*Var_max);
M.position = SortPosition(M.position);

x = [size_init, M.position, rule_set_initial];
M.cost = CostFunction(x,data);

% random number in [1,nVar]
b_rand = randperm(nVar,1);

% % trial vector T
% if (b_rand<=CR*nVar)
%     T = M;
% else
%     % the firt particle is used to initiate the vector T for now.
%     % Vector T will be recreated in the main loop of EPSO
%     T = particle(1);
% end

%% logging the best cost over the iterations
% array to hold the best cost function value at each iteration
best_costs = zeros(max_iteration, 1);

% setting PSO weights
if not(parameters.dynamic_weights)
    w = parameters.Kw1;
    c1 = parameters.Kc1;
    c2 = parameters.Kc2;
end

%% Main Loop of EPSO
for it=1:max_iteration
    fprintf('\nparEPSO_OP it #%d running\n\n',it);
    
    % updating weights
    if parameters.dynamic_weights
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
    end
    
    % making mutual vector and trial vector
    F = F_min + (F_max - F_min) * rand();  % scaling factor
    
    r = randperm(nPop, 3);             % random numbers(r1, r2, and r3)
    b_rand = randperm(nVar,1);         % random number in [1,nVar]
    
    % Matlab could not understand the variable M inside the paar-for loop.
    % To solve the problem, I initiated the the vector M as a temporary
    % solution.
    
    % making sure i is not the index in this par-for loop.
    i =0;
    dummy.position = particle(1).position;
    dummy.cost = particle(1).cost;
%     tic
    parfor z=1:nPop
        % update velocity
        particle(z).velocity = w*particle(z).velocity ...
            + c1*rand(Var_size).*(particle(z).best.position - particle(z).position) ...
            + c2*rand(Var_size).*(global_best.position - particle(z).position);
        
        % apply velocity limits
        particle(z).velocity = max(particle(z).velocity, Min_velocity);
        particle(z).velocity = min(particle(z).velocity, Max_velocity);
        
        % update position
        particle(z).position = particle(z).position + particle(z).velocity;
        
        % apply lower bound and upper ound limits
        particle(z).position = max(particle(z).position, Var_min);
        particle(z).position = min(particle(z).position, Var_max);
        
        % sorting particle(i).position for FLC parameters
        particle(z).position = SortPosition(particle(z).position);
        
        % evaluation
        x = [size_init, particle(z).position, rule_set_initial];
        [particle(z).cost, detail] = CostFunction(x,data);
            % particles_cost(i) = particle(i).cost;
        
        %% making mutual vector and trial vector
        % mutual vector M
        M = dummy;
        rnd = rand();
        M.position = Var_min + rnd*(Var_max - Var_min);
        
        % test this line to make sure that the max and min work over each element not the whole vector
        % apply lower bound and upper ound limits
        M.position = max(M.position, Var_min);
        M.position = min(M.position, Var_max);
        
        M.position = SortPosition(M.position);
        
        x = [size_init, M.position, rule_set_initial];
        [M.cost, M_detail] = CostFunction(x,data);
        
        M_flag = Constraints(x,data,M_detail);
        
        % random number in [1,nVar]
        b_rand = randperm(nVar,1);
        
        %         % trial vector T
        %         if (b_rand<=CR*nVar) && M_flag
        %             T = M;
        %         else
        %             % the firt particle is used to initiate the vector T for now.
        %             % Vector T will be recreated in the main loop of EPSO
        %             T = particle(z);
        %         end
        %         
        %         % improving position based on the trial vector
        %         if (T.cost<particle(z).cost)
        %             particle(z).position = T.position;
        %             particle(z).cost = T.cost;
        %         end
        
        % improving position based on the trial vector
        if (b_rand<=CR*nVar) && M_flag
            particle(z).position = M.position;
            particle(z).cost = M.cost;
        end
        
        % solution constraint check
        particle(z).flag = Constraints(x,data,detail);
        
        %         warning('particle flag only works because particle is indexed.')
        % update personal best
        if (particle(z).cost < particle(z).best.cost) && particle(z).flag
            particle(z).best.position = particle(z).position;
            particle(z).best.cost = particle(z).cost;
        end
    end
%     toc
      
    %     % update personal best
    %     parfor i=1:nPop
    %         if particle(i).cost<particle(i).best.cost
    %             particle(i).best.position = particle(i).position;
    %             particle(i).best.cost = particle(i).cost;
    %         end
    %     end

    % update global best
    tmp = [particle.best];
    [PB, GB_index] = sort([tmp.cost],'ascend');
        % global_best = tmp(GB_index(1));

    for counter=1:nPop
        if (PB(counter)<global_best.cost) && particle(GB_index(counter)).flag
            global_best = tmp(GB_index(counter));
            % calculate optimal rule_set for the new global optimum
            OR.stuck=1;
            try_it = 1;
            % optimal rule_set
            while (OR.stuck==1) && (try_it<parameters.try_it)
                disp('GB improved');
                OR = parEPSO_OR(problem, parameters, data, [size_init, global_best.position]);
                
                % updating global_best              
                if OR.best_solution.cost<global_best.cost
                    global_best.cost = OR.best_solution.cost;
                    % global_best position(FLC parameters stays the same, however
                    % the FLC rule_set gets updated).
                    % global_best.position = global_best.position;
                end
                   
                try_it = try_it + 1;
            end
            break 
        elseif (PB(counter)>=global_best.cost)
            break            
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
    blanked = full - length(size_init) - length(out.best_costs) -length(out.optimal_parameters) -length(out.optimal_rule_set) -7;
    new_file = [last_file;
                1111*ones(1,full);
                datetime_value,...
                out.stuck,data.w,parameters.max_iteration,...
                parameters.nPop,parameters.stuck_step,parameters.try_it,...
                out.best_costs',zeros(1,blanked),size_init,out.optimal_parameters,...
                out.optimal_rule_set;
                1111*ones(1,full)];
    xlswrite(filename,new_file,sheet,xlRange);
end

disp('parEPSO_OP ended.')
end

