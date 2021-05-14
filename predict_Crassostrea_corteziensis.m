function [prdData, info] = predict_Crassostrea_corteziensis(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
   if f_guzm >1 || f_guzm <0 || f_tL2 <0 || f_tL2 >1
    info = 0; prdData = {}; return;
   end  
  
   
  pars_T = [T_A T_H T_AH];
  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, pars_T);
  TC_tj = tempcorr(temp.tj, T_ref, pars_T);
  TC_tp = tempcorr(temp.tp, T_ref, pars_T);
  TC_am = tempcorr(temp.am, T_ref, pars_T);
  TC_Ri = tempcorr(temp.Ri, T_ref, pars_T);
  TC_tL = tempcorr(temp.tL, T_ref, pars_T);
  TC_tL2 = tempcorr(temp.tL2, T_ref, pars_T);
  TC_guzm = tempcorr(temp.F_T, T_ref, pars_T);
  
  % zero-variate data

  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve

  % hatch   
%   [U_H aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
%   a_h = aUL(2,1); 
%   aT_h = a_h/ TC_ah; % d, age at hatch at f and T
%   Lw_h = aUL(2,3)/ del_Me;           % cm, total length at hatch
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_M;                % cm, total length at birth at f
  Ww_b = L_b^3 *(1 + f * w);        % g, wet weight at birth
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth at f and T

  % metam
  L_j = L_m * l_j;                  % cm, structural length at metam
  Lw_j = L_j/ del_M;                % cm, total length at metam at f
  Ww_j = L_j^3 *(1 + f * w);        % g, wet weight at metam
  tT_j = (t_j - t_b)/ k_M/ TC_tj;   % d, time since birth at metam
  s_M = l_j/l_b;
  
  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, total length at puberty at f
  Ww_p = L_p^3 *(1 + f * w);        % g, wet weight at puberty 
  tT_p = (t_p - t_b)/ k_M/ TC_tp;   % d, time since birth at puberty at f and T

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb;U_Hj;  U_Hp]; % compose parameter vector at T
  RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  % pack to output
  prdData.ab = aT_b;
  prdData.tj = tT_j;
  prdData.tp = tT_p;
  prdData.am = aT_m;
  %prdData.LO = Lw_O;
  prdData.Lb = Lw_b;
  prdData.Lj = Lw_j;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Ri = RT_i;
  
  %% uni-variate data
  
%   % time-weight 
  t = tL;
  t(:,1) = t(:,1) + t_start;
  
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tL);
  kT_M = k_M * TC_tL; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_bj = L_b * exp(t((t(:,1) <= tT_j),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t((t(:,1) > tT_j),1) - tT_j)); % cm, expected length at time
  L = [L_bj; L_jm];
  ELw = L /del_M;
  
    % time-weight   
    t = tL2;
    t(:,1) = t(:,1) + t_start;
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tL2);
  kT_M = k_M * TC_tL2; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;        
  L_b = L_m * l_b;  L_j = L_m * l_j; L_i = L_m * l_i; tT_j = (t_j - t_b)/ kT_M;
  L_bj = L_b * exp(t((t(:,1) <= tT_j),1) * rT_j/ 3);
  L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t((t(:,1) > tT_j),1) - tT_j)); % cm, expected length at time
  L = [L_bj; L_jm];
  ELw2 = L /del_M;
  
  
  % Filtration

    L_guzm = (1/(f_guzm*w))^(1/3);            % cm, initial structural length
    EF_T = L_guzm^2  .* F_m .* f_guzm .* TC_guzm .* s_M./24;   %L/h, filtration rate           

  % Ingestion
    J_X = L_guzm^2 .* J_X_Am .* f_guzm .* TC_guzm .*s_M;   % mol/d, ingestion rate limited by salinity
    EJX_T = J_X .* w_X  /24 /d_X ;                  % g/h.g, ingestion rate dry
    
  % Oxygen consumption
    pars_pow = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp];
    p_ref = p_Am * L_m^2 * s_M;        % max assimilation power at max size
    pACSJGRD = p_ref * scaled_power_j(L_guzm, f_guzm, pars_pow, l_b, l_j, l_p);
    pADG = pACSJGRD(:, [1 7 5]);  
    O2M = (- n_M\n_O)';        % -, matrix that converts organic to mineral fluxes  O2M is prepared for post-multiplication eq. 4.35
    JO = pADG * eta_O';        %mol/d, organic fluxes
    JM = JO * O2M;             %mol/d, mineral fluxes
    EJO = -2*16 * JM(:,3) .* TC_guzm/24*1e3; % mgO2/h, oxygen consumption rate at temp T
    %EJN = 17 * JM(:,4) .* TC_guzm/24*1e6;   % ugNH4/h, nitrogen excretion rate at temp T
    
    
  % pack to output
  prdData.tL = ELw;
  prdData.tL2 = ELw2;
  prdData.F_T = EF_T;
  prdData.JX_T = EJX_T;
  prdData.JO_T= EJO;
  %prdData.JN_T= EJN;
  
  