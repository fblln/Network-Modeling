% FABIO ELLENA, LAB 1
clc
close all
clear all

len = 1000000; % number of jobs, if higher, the system is more similar to the theoric one
P_1 = 0.98; % percent rate of job 1 -> short
P_2 = 0.02; % percent rate of job 2 -> large
E_S_1 = 1; % size job 1
E_S_2 = 201; % size job 2

E_S = P_1 * E_S_1 + P_2 * E_S_2; % mean response time
E_S_second = P_1 * E_S_1 ^ 2 + P_2 * E_S_2 ^ 2; % mean second moment
rho_values = (0.1 : 0.1 : 0.9); % rho values

E_T_prevision = (rho_values ./ (1 - rho_values)) * (E_S_second) / (2 * E_S) + E_S; % prevision delay, calculated with pk + E_S
E_T_simulation = zeros(1, length(rho_values)); % simulated delay

job_sizes = rand(1, len); % size for each job that arrives
for i = 1:len
    if job_sizes(i) < P_1
        job_sizes(i) = E_S_1;
    else
        job_sizes(i) = E_S_2;
    end
end

counter = 0;
for counter = 1:length(rho_values)
    rho = rho_values(counter);
    lambda = rho / E_S; % calculate lambda, because rho = lambda * E_S
    
    arrival_interval = exprnd(1 / lambda, 1, len); % time passed until the next arrival
    arrival_time = cumsum(arrival_interval); % arrival time of each job
    service_time = zeros(1, len); % time at which each job is served
    delay = zeros(1, len); % time passed in the system by each job
    
    simulation_time = 0;
    for current_job = 1:len 
        % the idea is to artificially move the time until we find a job in
        % the system, then I register the output time for each job, and
        % later I use it to calculate the mean
        if simulation_time < arrival_time(current_job)
            simulation_time = arrival_time(current_job);
        end
        
        simulation_time = simulation_time + job_sizes(current_job); % job is served, time advances
        service_time(current_job) = simulation_time;
        delay(current_job) = service_time(current_job) - arrival_time(current_job);
    end
    
    E_T_simulation(counter) = mean(delay); % mean estimated delay
end

plot(rho_values, E_T_simulation, rho_values, E_T_prevision)
xlabel('\rho')
ylabel('E[T]')