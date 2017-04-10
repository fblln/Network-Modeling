% FABIO ELLENA, LAB 3
clc
close all
clear all

% The idea is to split the capacity of the server in two, in this way we
% create two 'virtual server' that serve only jobs with a specific size.
% short jobs are served by server1 and large jobs are served by server 2.
% The two server put together have the same capacity of the original one,
% but each job will be served at a lower speed.
% This is not important, because in this way jobs in each
% queue will have the same size, hence variability equal to zero, which
% meance low queuing time. This is because in this case short jobs are not
% stuck behind large ones, while in SJF this can happen.
% Each queue is an M/D/1, FCFS system, thus is quite easy to calculate the
% expected time in the system.
% The main problem is how to split the server capacity:
% The only hard constraint is that each queue should be stable, and I chose
% the values that meet this constraint when rho is high.
% Knowing the size distribution and the arrival rate of jobs, I can
% calculate the first derivative and find the optimal capacity split.

len = 1000000; % number of jobs, if higher, the system is more similar to the theoric one
P_1 = 0.98; % percent rate of job 1 -> short
P_2 = 0.02; % percent rate of job 2 -> large
E_S_1 = 1; % size job 1
E_S_2 = 201; % size job 2
E_S = P_1 * E_S_1 + P_2 * E_S_2; % mean response time (only service)

%This values are chosen to meet the constraint of rho_values_queue for each
%queue in the worst case : rho=0.9
% Also small changes can lend to an unstable system because the load is
% very high and we can't allow to have one of the two server in an unstable
% state
C_1 = 0.20; % capacity for job 1
C_2 = 0.80; % capacity for job 2

rho_values = (0.1 : 0.1 : 0.9); % rho values
lambda_values = rho_values ./E_S; % lambda values

E_S_Q1 = E_S_1 / C_1; %mean response time of queue 1
E_S_Q2 = E_S_2 / C_2; %mean response time queue 2

%rho should be always < 1, otherwise we have an unstable system and queuing
%time explodes
rho_values_Q1 = lambda_values * P_1 * E_S_Q1; %rho of queue 1
rho_values_Q2 = lambda_values * P_2 * E_S_Q2; %rho of queue 2

E_T_Q1 = rho_values_Q1 ./ ((1-rho_values_Q1)*2) * E_S_Q1 + E_S_Q1; % prevision delay in queue 1, calculated with pk + E_S_Q1
E_T_Q2 = rho_values_Q2 ./ ((1-rho_values_Q2)*2) * E_S_Q2 + E_S_Q2; % prevision delay in queue 2, calculated with pk + E_S_Q2

E_T_prevision = P_1 * E_T_Q1 + P_2 * E_T_Q2; % prevision delay
E_T_simulation = zeros(1, length(rho_values)); % simulated delay

job_sizes = rand(1, len); % size for each job that arrives

% assign job size to each job
for i = 1:len
    if job_sizes(i) < P_1
        job_sizes(i) = E_S_1;
    else
        job_sizes(i) = E_S_2;
    end
end

index_1 = find(job_sizes == E_S_1); % index of small jobs
index_2 = find(job_sizes == E_S_2); % index of large jobs
job_sizes_1 = job_sizes(index_1); % job of queue 1
job_sizes_2 = job_sizes(index_2); % job of queue 2

% start simulation for each value of rho
counter = 0;
for counter = 1:length(rho_values)
    rho = rho_values(counter); % pick a rho
    lambda = lambda_values(counter); % pick lambda
    
    arrival_interval = exprnd(1 / lambda, 1, len); % time passed until the next arrival
    arrival_time = cumsum(arrival_interval); % arrival time of each job
    
    arrival_time_1 = arrival_time(index_1); % arrival time of each job
    arrival_time_2 = arrival_time(index_2); % arrival time of each job
    
    service_time_1 = zeros(1, length(index_1)); % time at which each job is served
    service_time_2 = zeros(1, length(index_2)); % time at which each job is served
    
    delay_1 = zeros(1, length(arrival_time_1)); % time passed in the system by each job
    delay_2 = zeros(1, length(arrival_time_2)); % time passed in the system by each job
    
    
    % start simulation for system 1
    simulation_time = 0;
    
    % the idea is to artificially move the time until we find a job in
    % the system, then I register the output time for each job, and
    % later I use it to calculate di mean. I do this separately for the
    % queue 1 and the queue 2, which are separated. Then at the end I join
    % the results
    
    for current_job = 1:length(arrival_time_1)
        if simulation_time < arrival_time_1(current_job)
            simulation_time = arrival_time_1(current_job);
        end
        
        simulation_time = simulation_time + job_sizes_1(current_job)/ C_1; % job is served in Size/capacity because in this case, to each queue I allocate a percentage of the total capacity
        service_time_1(current_job) = simulation_time;
        delay_1(current_job) = service_time_1(current_job) - arrival_time_1(current_job); % time in the system
    end
    
    % start simulation for system 2
    simulation_time = 0;
    for current_job = 1:length(arrival_time_2) 
        if simulation_time < arrival_time_2(current_job)
            simulation_time = arrival_time_2(current_job);
        end
        
        simulation_time = simulation_time + job_sizes_2(current_job)/ C_2; % job is served in Size/capacity because in this case, to each queue I allocate a percentage of the total capacity
        service_time_2(current_job) = simulation_time;
        delay_2(current_job) = service_time_2(current_job) - arrival_time_2(current_job); % time in the system
    end
    
    delay = horzcat(delay_1,delay_2); % merge the delay arrays
    E_T_simulation(counter) = mean(delay); % mean estimated delay
end

plot(rho_values, E_T_simulation, rho_values, E_T_prevision)
xlabel('\rho')
ylabel('E[T]')