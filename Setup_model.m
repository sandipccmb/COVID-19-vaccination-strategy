clear all; 

state = 'India';

gps.age  = {'a1',  'a2', 'a3'};        % age groups                                
gps.risk = {'r0',  'r1'};              % risk groups
gps.vac  = {'v0',  'v1'};              % vaccination group
states1  = {'U','E','A','P','MS','R'}; % Variables
[i, s, d, lim] = get_addresses({states1, gps.vac, gps.risk,gps.age}, [], [], [], 0);
d = char(d);


% --- Include the auxiliaries ---------------------------------------------

n_agerisk = length(gps.risk)*length(gps.age);
n_age     = length(gps.age);
names     =     {'inc',      'mort'};      
lgths     = [n_agerisk,  n_agerisk];      
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;

s.infectious = [s.A,  s.P,  s.MS];                              
s.prevalent  = s.MS; 


% --- Incidence  --------------------------------------

ii = 1; inc_inds = {}; 
for ir = 1:length(gps.risk)
    for ia = 1:length(gps.age)
        inc_inds{ii} = intersect(intersect(s.(gps.risk{ir}), s.(gps.age{ia})),[s.MS]);
     ii = ii+1;
    end
end
[sel.inc,  agg.inc]  = get_aggregators(inc_inds,i.nstates);

% --- Counting deaths
mat = zeros(n_agerisk,i.nstates); row = 1; 
for ir = 1:length(gps.risk)
    for ia = 1:length(gps.age)
        mat(row, [intersect(intersect(s.MS, s.(gps.age{ia})), s.(gps.risk{ir}))]) = 1;
        row = row+1;
    end    
end
agg.mor = sparse(mat);  



% --- Input vector --------------------------------------------------------

fnames = {'incub_period','presym_period','infec_period','p_sympto','p_c','vmu','vcm'};
lengths= [1,1,1,1,1,6,9];
xi = []; lim = 0;
for ii = 1:length(fnames)
    inds = lim + [1:lengths(ii)];
    xi.(fnames{ii}) = inds;
    lim = inds(end);
end
bds = zeros(length(fnames),2);
bds(xi.incub_period,:)      = [3, 5];                                      % Uncertainity range of incubation period
bds(xi.presym_period,:)     = [0.5, 2];                                    % Uncertainity range of pre-symptomatic period
bds(xi.infec_period,:)      = [4, 7];                                      % Uncertainity range of infection period
bds(xi.p_sympto,:)          = [1/3, 2/3];                                  % Uncertainity range of proportion of symptomatic
bds(xi.p_c,:)               = [2/3, 1];                                    % Uncertainity range of relative infectiousness of asymptomatics
bds(xi.vmu,:)               = repmat([-1/4, 1/4],6,1);                     % variation in mortality rate
bds(xi.vcm,:)               = repmat([-1/4, 1/4],9,1);                     % variation in the element of contact matrix

prm.bounds                  = bds';


% --- Set up the parameters (CFR)among age groups (0-24y,25-60,60+) -----------------------------------------------
p.cfr          = [[0.001, 0.0145, 0.1092],2.5*[0.001, 0.0145, 0.1092]];    % Risk group-specific case fatality rates 

% --- Set up parameter values [vaccine effectiveness]
p.c1         = 0.6; %0.6;%0.9  [choose as required]                        % Relative susceptibility of the vaccinated groups
p.c2         = 0;                                                          % Relative infectiousness of the vaccinated groups compared to others 
p.c3         = 0; %0.6; 0.9    [choose as required]                        % Relative probability of developing symptomatic cases in the vaccinated groups

% Inter-age contact matrix 
prm.contact = [[1.3710,1.4317,0.0500];[2.5175,2.9007,0.0993];[0.2789,0.3350,0.0160]];% %Age group(0-24,25-60,60+) 
p.imm       = 1;

% --- Values that are overwritten by alloc_parameters
incub_period  = 4;                                                         % Average incubation period (days)
presym_period = 1;                                                         % Average pre-symptomatic period (days)
infec_period  = 5;                                                         % Average infectious period (days)

r.beta    = 1;
r.incub   = 1/incub_period;                                                % Rate of progression to infectiousness
p.sympto  = 2/3;                                                           % Proportion developing symptomatic infection
p.c       = 2/3;                                                           % Relative infectiousness of asymptomatic vs symptomatic infection
r.eta     = 1/presym_period;                                               % Rate from pre-symptomatic to symptomatic infection
r.gamma   = 1/infec_period;                                                % Rate of cure
r.mu      = p.cfr*r.gamma./(1-p.cfr);                                      % Mortality hazard

% --- Get population parameters -------------------------------------------
% Age-group of Population in each compartment [0-24y, 25-60y, >60]
prm.N     = [0.46, 0.445, 0.095]*1e6;  

% Proportion seropositive - later, do estimation based on final AR calculations
p.seropos = 0*ones(1,3);                        % keep zero here, adjust it in the initial condition
p.vacc1   = 0*ones(1,3);                        % Vaccination among non-risk group
p.vacc2   = 0*ones(1,3);                        % Vaccination among risk group

p.risk    = [0.028,0.143,0.43];                 %Proportion having comorbidity among different age groups[0-24y, 25-60y, >60]   

r.beta = 1;
M = make_model2(p, r, i, s, gps, prm);

inds = [s.E, s.infectious];

% Construct the new entries to E due to infection
tmp = repmat(M.lam,length(gps.risk)*length(gps.vac),1);
tmp((size(tmp,1)/2+1):end,:) = tmp((size(tmp,1)/2+1):end,:)*p.c1;

m = zeros(i.nstates);
m(s.E,:) = tmp;
% Account for the initial number in each compartment
mat = zeros(length(gps.age),length(gps.risk), length(gps.vac));
tmp = prm.N.*(1-p.seropos);
mat(:,1,1) = tmp.*(1-p.risk).*(1-p.vacc1);
mat(:,1,2) = tmp.*(1-p.risk).*p.vacc1;
mat(:,2,1) = tmp.*p.risk.*(1-p.vacc2);
mat(:,2,2) = tmp.*p.risk.*p.vacc2;
% Collapse them into vector form that lines up with elements of s.U
m(s.E,:) = m(s.E,:).*repmat(mat(:),1,size(m,2));               
F = m(inds,inds);
% Construct V
tmp = -(M.lin - diag(M.mortvec));
V   = tmp(inds,inds);
% Bring them together
r.beta = 1/max(eigs(V\F));

% % For a given R0 find beta value
% R0 = 2.5;
% r.beta1 = R0/max(eigs(V\F));


save Model_setup_India;
