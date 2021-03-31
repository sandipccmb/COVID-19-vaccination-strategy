% v2: Version to reflect updated model, i.e. with urbrur, age groups, and risk groups

function M = make_model2(p, r, i, s, gps, prm)


% -------------------------------------------------------------------------
% --- Get the linear rates ------------------------------------------------

m = zeros(i.nstates);
for iv = 1:length(gps.vac)
    vac = gps.vac{iv};
    for ia = 1:length(gps.age)
        age = gps.age{ia};
        for ir = 1:length(gps.risk)
            risk = gps.risk{ir};

               getaddr = @(st) i.(st).(vac).(risk).(age);
                
                U  = getaddr('U');                                         % Uninfected
                E  = getaddr('E');                                         % Exposed
                A  = getaddr('A');                                         % Asymptomatic
                P  = getaddr('P');                                         % Presymptomatic
                MS = getaddr('MS');                                        % Mild symptomatic
                R  = getaddr('R');                                         % Recovered
                
                % --- Development of symptoms (Mild or Severe)
                source  = P;
                destin  = MS;
                rate    = r.eta;
                m(destin, source) = m(destin, source) + rate;
                               
                % --- Recovering from disease
                sources = [A         MS];
                destin  = R;
                rate    = r.gamma;
                m(destin, sources) = m(destin, sources) + rate;              

        end
    end
end

for ia = 1:length(gps.age)
        age = gps.age{ia};
        for ir = 1:length(gps.risk)
            risk = gps.risk{ir};
            
                % --- Incubation
                source  = i.E.v0.(risk).(age);
                destins = [i.A.v0.(risk).(age),      i.P.v0.(risk).(age)];
                rates   = [(1-p.sympto),             p.sympto]*r.incub;
                m(destins, source) = m(destins, source) + rates';
        end
end

for ia = 1:length(gps.age)
        age = gps.age{ia};
        for ir = 1:length(gps.risk)
            risk = gps.risk{ir};
            
                % --- Incubation
                source  = i.E.v1.(risk).(age);
                destins = [i.A.v1.(risk).(age),   i.P.v1.(risk).(age)];
                rates   = [1-p.sympto*(1-p.c3),   p.sympto*(1-p.c3)]*r.incub;
   
                m(destins, source) = m(destins, source) + rates';
        end
end


M.lin = sparse(m - diag(sum(m,1)));


% -------------------------------------------------------------------------
% --- Get the nonlinear rates ---------------------------------------------

for ia = 1:length(gps.age)
    age = gps.age{ia};
      
        m = zeros(i.nstates);
        for iv = 1:length(gps.vac)
            vac = gps.vac{iv};
            for ir = 1:length(gps.risk)
                risk = gps.risk{ir};
                
                getaddr = @(st) i.(st).(vac).(risk).(age);
                U  = getaddr('U');                                         % Uninfected
                E  = getaddr('E');                                         % Exposed
                R  = getaddr('R');                                         % Recovered
                
                m(E, U) = 1;
                m(E, R) = 1-p.imm;
            end
        end
        
        m(:,s.v1) = m(:,s.v1)*(1-p.c1);
        M.nlin.(age) = sparse(m - diag(sum(m,1)));
%    end
end

% -------------------------------------------------------------------------
% --- Getting force-of-infection ------------------------------------------
% p.c1: Relative susceptibility
% p.c2: Relative infectiousness
% p.c: Relative infectiousness of asymptomatic vs symptomatic infection

% First, construct basic building block of repeating contact matrix, that can be weighted as necessary
templ = zeros(3,i.nstates);
templ(:,s.infectious) = repmat(prm.contact,1,length(s.infectious)/length(gps.age));

% Adjust for asymptomatics being less infectious
templ(:,s.A) = templ(:,s.A)*p.c;

% % Adjust for vaccinated people being less infectious
% templ(:,s.v1) = templ(:,s.v1)*(1-p.c2);

% Correct for population numbers
tmp2  = prm.N(:)';
den   = repmat(tmp2, size(templ,1), length(s.infectious)/length(gps.age)); 
templ(:,s.infectious) = templ(:,s.infectious)./den;

m = r.beta*templ;
M.lam = sparse(m);

% --- Get the mortality rates
m = zeros(1,i.nstates);
ii = 1;
for ir = 1:length(gps.risk)
    risk = gps.risk{ir};
    for ia = 1:length(gps.age)
        age = gps.age{ia};
        inds = intersect(s.MS,intersect(s.(risk),s.(age)));
        m(inds) = r.mu(ii); ii = ii+1;
    end
end
M.mortvec = m';
