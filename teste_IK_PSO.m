%% EXPERIMENTO DE CINEMÁTICA INVERSA (MONTE CARLO COM CRONÔMETRO)
% Comparação entre PSO Clássico e M-FRPSO para o manipulador Pioneer 7DOF
clc;
clear;
close all;

%% 1. PARÂMETROS DO ROBÔ E LIMITES DE JUNTAS (Pioneer 7DOF)
% Limites de junta (radianos)
qmin = deg2rad([-150, -88, -150, -92.5, -150, -105.5, -102.5])';
qmax = deg2rad([ 150,  88,  150,  92.5,  150,  105.5,  102.5])';

%% 2. PARÂMETROS GERAIS DA SIMULAÇÃO
N_ENSAIOS = 1000;    % Número de inicializações estocásticas
MAX_ITER  = 1000;    % Limite máximo de iterações por ensaio
TOL_ERRO  = 0.001;   % Tolerância para considerar convergência

% Pesos da função de custo
wd = 1.0; % Peso para distância de posição
wo = 0.1; % Peso para distância de orientação

%% 3. ESTRUTURAS DE ARMAZENAMENTO DE RESULTADOS
% Guardamos se convergiu, número de iterações e o tempo (em segundos)
res.PSO_Classico.TarefaA.convergidos = false(N_ENSAIOS, 1);
res.PSO_Classico.TarefaA.iteracoes   = zeros(N_ENSAIOS, 1);
res.PSO_Classico.TarefaA.tempos      = zeros(N_ENSAIOS, 1);

res.PSO_Classico.TarefaB.convergidos = false(N_ENSAIOS, 1);
res.PSO_Classico.TarefaB.iteracoes   = zeros(N_ENSAIOS, 1);
res.PSO_Classico.TarefaB.tempos      = zeros(N_ENSAIOS, 1);

res.M_FRPSO.TarefaA.convergidos      = false(N_ENSAIOS, 1);
res.M_FRPSO.TarefaA.iteracoes        = zeros(N_ENSAIOS, 1);
res.M_FRPSO.TarefaA.tempos           = zeros(N_ENSAIOS, 1);

res.M_FRPSO.TarefaB.convergidos      = false(N_ENSAIOS, 1);
res.M_FRPSO.TarefaB.iteracoes        = zeros(N_ENSAIOS, 1);
res.M_FRPSO.TarefaB.tempos           = zeros(N_ENSAIOS, 1);

%% 4. EXECUÇÃO DO LAÇO DE MONTE CARLO (1000 ENSAIOS)
fprintf('==================================================\n');
fprintf('Iniciando Experimentos de Monte Carlo (%d ensaios)\n', N_ENSAIOS);
fprintf('==================================================\n');

% INÍCIO DA CONTAGEM DE TEMPO TOTAL
tempo_total_inicio = tic;

for ensaio = 1:N_ENSAIOS
    if mod(ensaio, 100) == 0 || ensaio == 1
        fprintf('Executando Ensaio %d de %d...\n', ensaio, N_ENSAIOS);
    end
    
    % Sorteia uma pose de destino alcançável a partir de uma configuração válida
    q_alvo = qmin + rand(7, 1) .* (qmax - qmin);
    T_alvo = fkin_pioneer(q_alvo);
    p_alvo = T_alvo(1:3, 4);
    o_alvo = eulerZYX(T_alvo(1:3, 1:3));
    
    % --- TAREFA A: APENAS POSIÇÃO ---
    % 1. PSO Clássico
    t_start = tic;
    [conv, iters] = rodar_PSO_Classico(p_alvo, o_alvo, qmin, qmax, MAX_ITER, TOL_ERRO, 1.0, 0.0);
    res.PSO_Classico.TarefaA.tempos(ensaio) = toc(t_start);
    res.PSO_Classico.TarefaA.convergidos(ensaio) = conv;
    res.PSO_Classico.TarefaA.iteracoes(ensaio)   = iters;
    
    % 2. M-FRPSO
    t_start = tic;
    [conv, iters] = rodar_M_FRPSO(p_alvo, o_alvo, qmin, qmax, MAX_ITER, TOL_ERRO, 1.0, 0.0);
    res.M_FRPSO.TarefaA.tempos(ensaio) = toc(t_start);
    res.M_FRPSO.TarefaA.convergidos(ensaio)      = conv;
    res.M_FRPSO.TarefaA.iteracoes(ensaio)        = iters;
    
    % --- TAREFA B: POSE COMPLETA ---
    % 1. PSO Clássico
    t_start = tic;
    [conv, iters] = rodar_PSO_Classico(p_alvo, o_alvo, qmin, qmax, MAX_ITER, TOL_ERRO, wd, wo);
    res.PSO_Classico.TarefaB.tempos(ensaio) = toc(t_start);
    res.PSO_Classico.TarefaB.convergidos(ensaio) = conv;
    res.PSO_Classico.TarefaB.iteracoes(ensaio)   = iters;
    
    % 2. M-FRPSO
    t_start = tic;
    [conv, iters] = rodar_M_FRPSO(p_alvo, o_alvo, qmin, qmax, MAX_ITER, TOL_ERRO, wd, wo);
    res.M_FRPSO.TarefaB.tempos(ensaio) = toc(t_start);
    res.M_FRPSO.TarefaB.convergidos(ensaio)      = conv;
    res.M_FRPSO.TarefaB.iteracoes(ensaio)        = iters;
end

% FIM DA CONTAGEM DE TEMPO TOTAL
tempo_total_decorrido = toc(tempo_total_inicio);

%% 5. PROCESSAMENTO E EXIBIÇÃO DOS RESULTADOS
fprintf('\n=================================================================================\n');
fprintf('                             RESULTADOS DOS EXPERIMENTOS                         \n');
fprintf('=================================================================================\n');

imprimir_metricas('PSO Clássico - Tarefa A', res.PSO_Classico.TarefaA);
imprimir_metricas('M-FRPSO      - Tarefa A', res.M_FRPSO.TarefaA);
fprintf('---------------------------------------------------------------------------------\n');
imprimir_metricas('PSO Clássico - Tarefa B', res.PSO_Classico.TarefaB);
imprimir_metricas('M-FRPSO      - Tarefa B', res.M_FRPSO.TarefaB);
fprintf('=================================================================================\n');
fprintf('Tempo total gasto para os %d ensaios: %.2f segundos (%.2f minutos)\n', ...
        N_ENSAIOS, tempo_total_decorrido, tempo_total_decorrido/60);
fprintf('=================================================================================\n');


%% =======================================================================
%% ALGORITMO 1: M-FRPSO (MODIFIED FULLY RESAMPLED PSO)
%% =======================================================================
function [convergiu, iter] = rodar_M_FRPSO(p_alvo, o_alvo, qmin, qmax, max_iter, tol, wd, wo)
    P = 100;     % Tamanho da população
    N = 5;       % Partículas de Elite
    tau = 0.5;   % Parâmetro de escala da dispersão
    u = (P - N)/N;
    
    Pop = qmin + rand(7, P) .* (qmax - qmin);
    J = zeros(P, 1);
    
    convergiu = false;
    for iter = 1:max_iter
        for i = 1:P
            T = fkin_pioneer(Pop(:, i));
            p_atual = T(1:3, 4);
            o_atual = eulerZYX(T(1:3, 1:3));
            
            dp = norm(p_alvo - p_atual);
            do = norm(o_alvo - o_atual);
            J(i) = wd * dp + wo * do;
        end
        
        [val, idx] = sort(J);
        Pop = Pop(:, idx);
        
        if val(1) < tol
            convergiu = true;
            return;
        end
        
        Popbests = Pop(:, 1:N);
        sigma = val(1) / tau;
        
        for i = 1:N
            Popbest = Popbests(:, i);
            Pop(:, N + 1 + u*(i-1) : N + u*i) = sigma * randn(7, u) + Popbest;
        end
        
        Pop = max(min(Pop, qmax), qmin);
    end
end

%% =======================================================================
%% ALGORITMO 2: PSO CLÁSSICO (INERTIA WEIGHT PSO)
%% =======================================================================
function [convergiu, iter] = rodar_PSO_Classico(p_alvo, o_alvo, qmin, qmax, max_iter, tol, wd, wo)
    P = 100;        % Tamanho da população
    w = 0.7;        % Inércia
    c1 = 1.5;       % Cognitivo
    c2 = 1.5;       % Social
    
    Pop = qmin + rand(7, P) .* (qmax - qmin);
    V   = zeros(7, P);
    
    Pbest = Pop;
    J_pbest = inf(P, 1);
    
    Gbest = Pop(:, 1);
    J_gbest = inf;
    
    convergiu = false;
    for iter = 1:max_iter
        for i = 1:P
            T = fkin_pioneer(Pop(:, i));
            p_atual = T(1:3, 4);
            o_atual = eulerZYX(T(1:3, 1:3));
            
            dp = norm(p_alvo - p_atual);
            do = norm(o_alvo - o_atual);
            custo = wd * dp + wo * do;
            
            if custo < J_pbest(i)
                J_pbest(i) = custo;
                Pbest(:, i) = Pop(:, i);
            end
            
            if custo < J_gbest
                J_gbest = custo;
                Gbest = Pop(:, i);
            end
        end
        
        if J_gbest < tol
            convergiu = true;
            return;
        end
        
        r1 = rand(7, P);
        r2 = rand(7, P);
        V = w * V + c1 * r1 .* (Pbest - Pop) + c2 * r2 .* (Gbest - Pop);
        Pop = Pop + V;
        
        Pop = max(min(Pop, qmax), qmin);
    end
end

%% =======================================================================
%% FUNÇÕES AUXILIARES DE CINEMÁTICA E UTILITÁRIOS
%% =======================================================================

function T7 = fkin_pioneer(q)
    a1 = 0;     alpha1 =  pi/2; d1 = 0.125;  theta1 =  pi/2 + q(1);
    a2 = 0;     alpha2 = -pi/2; d2 = 0;      theta2 = q(2);
    a3 = 0;     alpha3 =  pi/2; d3 = 0.150;  theta3 = q(3);
    a4 = 0;     alpha4 = -pi/2; d4 = 0;      theta4 = q(4);
    a5 = 0;     alpha5 =  pi/2; d5 = 0.145;  theta5 = q(5);
    a6 = 0.075; alpha6 =  pi/2; d6 = 0;      theta6 = q(6) - pi/2;
    a7 = 0;     alpha7 =  pi/2; d7 = 0;      theta7 = q(7);

    A1 = Matriz_A(theta1, a1, alpha1, d1);
    A2 = Matriz_A(theta2, a2, alpha2, d2);
    A3 = Matriz_A(theta3, a3, alpha3, d3);
    A4 = Matriz_A(theta4, a4, alpha4, d4);
    A5 = Matriz_A(theta5, a5, alpha5, d5);
    A6 = Matriz_A(theta6, a6, alpha6, d6);
    A7 = Matriz_A(theta7, a7, alpha7, d7);

    T7 = A1 * A2 * A3 * A4 * A5 * A6 * A7;
end

function A = Matriz_A(theta, a, alpha, d)
    A = [ cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
          sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
          0,           sin(alpha),             cos(alpha),             d;
          0,           0,                      0,                      1 ];
end

function o = eulerZYX(R)
    if abs(R(3,1)) < 1
        pitch = -asin(R(3,1));
        roll  = atan2(R(3,2), R(3,3));
        yaw   = atan2(R(2,1), R(1,1));
    else
        pitch = sign(-R(3,1)) * pi/2;
        roll  = 0;
        yaw   = atan2(-R(1,2), R(2,2));
    end
    o = [yaw, pitch, roll];
end

function imprimir_metricas(nome, dados)
    n_conv = sum(dados.convergidos);
    tc = (n_conv / length(dados.convergidos)) * 100;
    
    if n_conv > 0
        mi = mean(dados.iteracoes(dados.convergidos));
        tempo_medio_ms = mean(dados.tempos(dados.convergidos)) * 1000; % Em milissegundos
    else
        mi = NaN;
        tempo_medio_ms = NaN;
    end
    
    fprintf('%-25s -> tc: %6.2f%% | mi: %6.2f it | t_médio: %6.2f ms\n', ...
            nome, tc, mi, tempo_medio_ms);
end