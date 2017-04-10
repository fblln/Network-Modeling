% FABIO ELLENA, LAB 2
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
lambda_values = rho_values/E_S; % lambda values

%every job is put in a virtual queue of jobs with the same size, and queues
%with jobs of smaller size have higher priority
E_T_Q1 = rho_values.*(E_S_second/(2*E_S))./(1-P_1*lambda_values*E_S_1);
E_T_Q2 = rho_values.*(E_S_second/(2*E_S))./(1-P_1*lambda_values*1-P_2*lambda_values*E_S_2)./(1-P_1*lambda_values*E_S_1);
E_T_prevision = P_1 * E_T_Q1 + P_2 * E_T_Q2 + E_S; % prevision delay, calculated with pk + E_S

E_T_simulation = zeros(1, length(rho_values)); % simulated delay

job_sizes = rand(1, len); % size for each job that arrives
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

counter = 0;
for counter = 1:length(rho_values)
    rho = rho_values(counter); % pick a rho
    lambda = rho / E_S; % pick lambda
    
    arrival_interval = exprnd(1 / lambda, 1, len); % time passed until the next arrival
    arrival_time = cumsum(arrival_interval); % arrival time of each job
    
    arrival_time_1 = arrival_time(index_1); % arrival time of each job
    arrival_time_2 = arrival_time(index_2); % arrival time of each job
    
    service_time_1 = zeros(1, length(index_1)); % time at which each job is served
    service_time_2 = zeros(1, length(index_2)); % time at which each job is served
    
    delay_1 = zeros(1, length(arrival_time_1)); % time passed in the system by each job
    delay_2 = zeros(1, length(arrival_time_2)); % time passed in the system by each job
    
    simulation_time = 0;
    current_1 = 1; %current job to process
    current_2 = 1; %current job to process
    len_1 = length(arrival_time_1); % number of short jobs
    len_2 = length(arrival_time_2); % number of large jobs
    
    %start simulation. Here I divided jobs into two virtual queues to
    %improve the processing time, which is O(N). But there is only one real queue and the SJF policy
    %is implemented
    
    %I have a counter for each virtual queue and a unique simulation_time.
    %the simulation time start from 0 and when there are no job that
    %arrived in the past, I pick the first one that will arrive in the
    %future and I fast forward the time. Otherwise (jobs in the past to
    %serve), I pick the earliest short job in queue and I serve it. If there
    %are not anymore short jobs, I pick a large one.
    
    while current_1 <= len_1 || current_2 <= len_2
        % try to execute a short job that is in the queue
        if current_1 <= len_1 && arrival_time_1(current_1)<= simulation_time
            simulation_time = simulation_time + job_sizes_1(current_1); % job is served
            service_time_1(current_1) = simulation_time;
            delay_1(current_1) = service_time_1(current_1) - arrival_time_1(current_1); % time in the system
            current_1 = current_1 + 1; % increment short job pointer
            % try to execute a large job that is in the queue
        elseif current_2 <= len_2 && arrival_time_2(current_2)<= simulation_time
            simulation_time = simulation_time + job_sizes_2(current_2); % job is served
            service_time_2(current_2) = simulation_time;
            delay_2(current_2) = service_time_2(current_2) - arrival_time_2(current_2); % time in the system
            current_2 = current_2 + 1; % increment short job pointer
        else
            % move forward the time, because there is no job in the queue
            %if there are not anymore short jobs, pick the first large job
            if current_1 > len_1
                simulation_time = arrival_time_2(current_2) + job_sizes_2(current_2); % job is served
                service_time_2(current_2) = simulation_time;
                delay_2(current_2) = service_time_2(current_2) - arrival_time_2(current_2); % time in the system
                current_2 = current_2 + 1; % increment large job pointer
                %if there are not anymore large jobs, pick the first small job
            elseif current_2 > len_2
                simulation_time = arrival_time_1(current_1) + job_sizes_1(current_1); % job is served
                service_time_1(current_1) = simulation_time;
                delay_1(current_1) = service_time_1(current_1) - arrival_time_1(current_1); % time in the system
                current_1 = current_1 + 1; % increment short job pointer
            else
                % There are short and large jobs in the future, find the
                % earliest one
                if arrival_time_1(current_1)<=arrival_time_2(current_2)
                    simulation_time = arrival_time_1(current_1) + job_sizes_1(current_1); % job is served
                    service_time_1(current_1) = simulation_time;
                    delay_1(current_1) = service_time_1(current_1) - arrival_time_1(current_1); % time in the system
                    current_1 = current_1 + 1; % increment short job pointer
                else
                    simulation_time = arrival_time_2(current_2) + job_sizes_2(current_2); % job is served
                    service_time_2(current_2) = simulation_time;
                    delay_2(current_2) = service_time_2(current_2) - arrival_time_2(current_2); % time in the system
                    current_2 = current_2 + 1; % increment short job pointer
                end
            end
        end
    end
    
    delay = horzcat(delay_1,delay_2); % merge the delay arrays
    E_T_simulation(counter) = mean(delay); % mean estimated delay
end

plot(rho_values, E_T_simulation, rho_values, E_T_prevision)
xlabel('\rho')
ylabel('E[T]')