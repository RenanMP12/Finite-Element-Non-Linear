clear all; close; clc;
%% step 1
n = 0; % step
t = 0; % time

%% step 2 - set all nodal displacements to zero
nel = 5;                    % number of elements
nnodes = nel*4 - (nel - 1); % total number of nodes
u = zeros(nnodes,1);        % global displacement vector     
v_half = zeros(nnodes,1);        % global velocity vector     

%% step 3 - initial stress BC
sigma_0 = zeros(nnodes,1);  % global initial stress  
BCs = [1 0];              % setting boundary condition [node | value]
alldof = linspace(1, nnodes, nnodes);
freedof = alldof;
for i = 1 : size(BCs,1)
    freedof(BCs(i,1)) = BCs(i,2);
end

%% step 4 - compute mass matrix
E = 70e9;                     % Youngs modulus
t_max = 0.01;                 % max time interval
alpha = 0.9;                  % value of alpha to get delta_t
b = 9.81;                     % gravity/body weight in m/s2
rho = 2700;                   % density in kg/m3
area = 300e-6;                % initial cross sectional area in m2
length = 200e-3;              % initial length of the rod in m
he = length/nel;              % elemental length
mass_e = zeros(4,4);          % elemental mass matrix
crit_t = 1e-6;                % critical time interval
a_n = zeros(nnodes, 1);       % acceleration matrix pre allocation

for i = 1 : 4 % Integration 4 Gauss points
    x_gauss = [-0.861136, -0.339981, 0.339981, 0.861136];
    w_gauss = [0.347855, 0.652145, 0.652145, 0.347855];
    x = x_gauss(i); w = w_gauss(i);
    
    % shape function and derivative
    N = 1/16*[-9*x^3 + 9*x^2 + x - 1, 27*x^3 - 9*x^2 - 27*x + 9, -27*x^3 - 9*x^2 + 27*x + 9, 9*x^3 + 9*x^2 - x - 1];
    mass_e =  mass_e + rho*area*N'*N*w*he/2;
end

% diagonalize mass matrix
mass_e = sum(mass_e(1,:))*[1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] +...
         sum(mass_e(2,:))*[0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0] +...
         sum(mass_e(3,:))*[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0] +...
         sum(mass_e(4,:))*[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1];

% assembly global mass matrix
mass_global = zeros(nnodes, nnodes);
for i = 1 : nel
    loc = (i - 1)*3 + [1 2 3 4];
    mass_global(loc, loc) = mass_global(loc, loc) + mass_e;
end

%% step 5 - compute force (stress update) subroutine

t_vector = zeros(1e4,1);
u_vector = zeros(1e4,1);
v_vector = zeros(1e4,1);
a_vector = zeros(1e4,1);
cond = false;
cycle_number = 1;

while cond == false
    % first sub-step
    t_vector(cycle_number) = t;
    u_vector(cycle_number) = u(end);
    v_vector(cycle_number) = v_half(end);
    a_vector(cycle_number) = a_n(end);

    local_displacement = zeros(nel, 4); % elemental displacement at step n
    local_velocity = zeros(nel, 4);     % elemental velocity at step n_half
    
    % second sub-step: Loop over elements e = 1 to ne
    
    % "gather" element nodal displacements and velocities
    for i = 1 : nel
        loc = (i - 1)*3 + [1 2 3 4];
        local_displacement(i, :) = u(loc);
        local_velocity(i, :) = v_half(loc);
    end
    
    % pre-allocating internal and external elemental force vectors
    f_int_n = zeros(4,1, nel);
    f_ext_n = zeros(4,1, nel);
    
    % force at the right tip
%     tau = 1000e3; label = 'Instataneous load'; % value given in N
    tau = 100000e3*t; label = 'Linear load'; % value given in N
    stiffness_e = zeros(4,4,nel); % elemental stiffness matrix
    % third sub-step: loop over Gauss points xi_q for q = 1 : nq
    for j = 1:nel
        % compute B
        for i = 1 : 3 % Integration 3 Gauss points
            x_gauss = [0.0000000000,-0.7745966692, 0.7745966692];
            w_gauss = [0.8888888889, 0.5555555555, 0.5555555555];
            x = x_gauss(i); w = w_gauss(i);
            
            % shape function and derivative
            N = 1/16*[-9*x^3 + 9*x^2 + x - 1, 27*x^3 - 9*x^2 - 27*x + 9, -27*x^3 - 9*x^2 + 27*x + 9, 9*x^3 + 9*x^2 - x - 1];
            B = [-1.6875*x^2 + 1.125*x + 0.0625, 5.0625*x^2 - 1.125*x - 1.6875, -5.0625*x^2 - 1.125*x + 1.6875, 1.6875*x^2 + 1.125*x - 0.0625];
            stiffness_e(:,:,j) =  stiffness_e(:,:,j) + E*area*B'*B*w*2/he;
            f_ext_n(:,:,j) = f_ext_n(:,:,j) + w*rho*N'*area*b*he/2;
        end
        if j == nel
            f_ext_n(:,:,j) = f_ext_n(:,:,j) + [0;0;0;tau];
        end
        f_int_n(:,1,j) = stiffness_e(:,:,j)*local_displacement(j,:)';
    end
    
    % forth sub-step: compute elemental force
    f_n = f_ext_n - f_int_n;
    
    % assembly global force vector
    f_n_global = zeros(nnodes, 1);
    for i = 1 : nel
        loc = (i - 1)*3 + [1 2 3 4];
        f_n_global(loc, 1) = f_n_global(loc, 1) + f_n(:,:,i);
    end
    
    % find critical time interval
    delta_t = alpha*crit_t;
    
    %% step 6 - compute acceleration
    
    a_n(logical(freedof)) = mass_global(logical(freedof),logical(freedof))\f_n_global(logical(freedof));
    
    %% step 7 & step 8 - nodal velocities at half steps & enforce the velocity BCs
    if n == 0
        v_half = v_half + a_n*delta_t/2;
    else
        v_half = v_half + a_n*delta_t;
    end
    
    %% step 9 & step 10 - update nodal displacements & enforce the displacement BCs
    u = u + v_half*delta_t;  
    
    %% step 11 - update time step n and time t
    t = t + delta_t;
    n = n + 1;
    
    %% step 12 - if the simulation time is not done ---> go to step 5
    cycle_number = cycle_number + 1;
    if t > t_max 
        cond = true;
    end
end

%% plot
% tiledlayout(3,1)

% Top plot
% ax1 = nexttile;
% plot(ax1, t_vector, u_vector)
figure(1)
hold on
plot(t_vector, u_vector)
% grid(ax1,'on')
grid on
ylabel('Displacement (m)','FontSize', 20)
set(gca,'FontSize',20)
legend('Linear load','Instataneous load')

% ax2 = nexttile;
% plot(ax2, t_vector, v_vector)
% ylabel('Velocity (m/s)','FontSize', 20)
% set(gca,'FontSize',20)
% 
% ax3 = nexttile;
% plot(ax3, t_vector, a_vector)
xlabel('Time (s)','FontSize', 20)
% ylabel('Acceleration (m/s^{2})','FontSize', 20)
% set(gca,'FontSize',20)
