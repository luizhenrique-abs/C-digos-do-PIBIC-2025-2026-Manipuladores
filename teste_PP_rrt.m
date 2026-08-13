%% EXPERIMENTO DE PLANEJAMENTO DE CAMINHO (MONTE CARLO)
% Comparação entre RRT Clássica e RRT Polarizada (DRRT) para o Pioneer 7DOF
clc;
clear;
close all;

%% 1. PARÂMETROS DO ROBÔ E LIMITES DE JUNTAS (Pioneer 7DOF)
qmin = deg2rad([-150, -88, -150, -92.5, -150, -105.5, -102.5])';
qmax = deg2rad([ 150,  88,  150,  92.5,  150,  105.5,  102.5])';

%% 2. PARÂMETROS ATUALIZADOS (Linhas 10 a 16 do seu script)
N_ENSAIOS  = 100;      % 100 cenários de planejamento
MAX_NODES  = 10000;    % Aumentado para dar mais chance à RRT Clássica
PASSO_D    = 0.15;     % Passo de expansão aumentado (0.15 rad)
GOAL_TOL   = 0.20;     % Tolerância para considerar meta alcançada
DELTA_MIN  = 0.08;     % Distância mínima para descarte de nós (DRRT)
PROB_POLAR = 0.25;     % Polarização de 25% na RRT Polarizada

%% 3. DEFINIÇÃO DE OBSTÁCULOS FISICOS NO ESPAÇO DE TRABALHO
% Obstáculo 1: Caixa/Mesa frontal
obs1.tipo = 'box'; obs1.bounds = [-0.4, 0.4; 0.2, 0.6; -0.2, 0.3]; 
% Obstáculo 2: Esfera simulando restrição/passagem estreita
obs2.tipo = 'sphere'; obs2.center = [0.0; 0.35; 0.25]; obs2.radius = 0.15;

obstaculos = {obs1, obs2};

%% 4. ESTRUTURAS DE ARMAZENAMENTO DE RESULTADOS
res.RRT_Classica.convergidos = false(N_ENSAIOS, 1);
res.RRT_Classica.tempos      = zeros(N_ENSAIOS, 1);
res.RRT_Classica.nos_totais  = zeros(N_ENSAIOS, 1);
res.RRT_Classica.comprimento = zeros(N_ENSAIOS, 1);

res.RRT_Polarizada.convergidos = false(N_ENSAIOS, 1);
res.RRT_Polarizada.tempos      = zeros(N_ENSAIOS, 1);
res.RRT_Polarizada.nos_totais  = zeros(N_ENSAIOS, 1);
res.RRT_Polarizada.comprimento = zeros(N_ENSAIOS, 1);

%% 5. EXECUÇÃO DOS EXPERIMENTOS DE MONTE CARLO
fprintf('======================================================================\n');
fprintf('Iniciando Experimentos de Planejamento de Caminho (%d ensaios)\n', N_ENSAIOS);
fprintf('======================================================================\n');

tempo_total_inicio = tic;

for ensaio = 1:N_ENSAIOS
    if mod(ensaio, 10) == 0 || ensaio == 1
        fprintf('Executando Ensaio de Planejamento %d de %d...\n', ensaio, N_ENSAIOS);
    end
    
    % Encontra um par (q_init, q_goal) livre de colisões
    [q_init, q_goal] = gerar_par_valido(qmin, qmax, obstaculos);
    
    % --- 1. RRT CLÁSSICA ---
    t_start = tic;
    [conv, nos, comp] = rodar_RRT_Classica(q_init, q_goal, qmin, qmax, ...
                                           obstaculos, MAX_NODES, PASSO_D, GOAL_TOL);
    res.RRT_Classica.tempos(ensaio)      = toc(t_start);
    res.RRT_Classica.convergidos(ensaio) = conv;
    res.RRT_Classica.nos_totais(ensaio)  = nos;
    res.RRT_Classica.comprimento(ensaio) = comp;
    
    % --- 2. RRT POLARIZADA (DRRT) ---
    t_start = tic;
    [conv, nos, comp] = rodar_RRT_Polarizada(q_init, q_goal, qmin, qmax, ...
                                             obstaculos, MAX_NODES, PASSO_D, ...
                                             GOAL_TOL, DELTA_MIN, PROB_POLAR);
    res.RRT_Polarizada.tempos(ensaio)      = toc(t_start);
    res.RRT_Polarizada.convergidos(ensaio) = conv;
    res.RRT_Polarizada.nos_totais(ensaio)  = nos;
    res.RRT_Polarizada.comprimento(ensaio) = comp;
end

tempo_total_decorrido = toc(tempo_total_inicio);

%% 6. EXIBIÇÃO DAS MÉTRICAS GLOBAIS E ESTRUTURAIS
fprintf('\n=============================================================================================\n');
fprintf('                             RESULTADOS DE PLANEJAMENTO DE CAMINHO                           \n');
fprintf('=============================================================================================\n');

imprimir_metricas_pp('RRT Clássica   ', res.RRT_Classica);
imprimir_metricas_pp('RRT Polarizada ', res.RRT_Polarizada);

fprintf('=============================================================================================\n');
fprintf('Tempo total dos %d ensaios de PP: %.2f segundos (%.2f minutos)\n', ...
        N_ENSAIOS, tempo_total_decorrido, tempo_total_decorrido/60);
fprintf('=============================================================================================\n');


%% =======================================================================
%% ALGORITMO 1: RRT CLÁSSICA (COM 5% DE BIAS PADRÃO)
function [convergiu, num_nos, comprimento] = rodar_RRT_Classica(q_init, q_goal, qmin, qmax, obstaculos, max_nodes, d, goal_tol)
    Arvore = [q_init; 0];
    convergiu = false;
    comprimento = NaN;
    
    for k = 1:max_nodes
        % Bias minimalista clássico de 5% apenas para permitir exploração em 7D
        if rand < 0.05
            q_rand = q_goal;
        else
            q_rand = qmin + rand(7, 1) .* (qmax - qmin);
        end
        
        dists = sqrt(sum((Arvore(1:7, :) - q_rand).^2, 1));
        [~, idx_near] = min(dists);
        q_near = Arvore(1:7, idx_near);
        
        dir = q_rand - q_near;
        norm_dir = norm(dir);
        if norm_dir == 0, continue; end
        q_new = q_near + d * (dir / norm_dir);
        q_new = max(min(q_new, qmax), qmin);
        
        if ~checa_colisao_segmento(q_near, q_new, obstaculos)
            Arvore = [Arvore, [q_new; idx_near]];
            if norm(q_new - q_goal) < goal_tol
                convergiu = true;
                break;
            end
        end
    end
    
    num_nos = size(Arvore, 2);
    if convergiu
        comprimento = calcula_comprimento_caminho(Arvore);
    end
end

%% =======================================================================
%% ALGORITMO 2: RRT POLARIZADA (DRRT)
function [convergiu, num_nos, comprimento] = rodar_RRT_Polarizada(q_init, q_goal, qmin, qmax, obstaculos, max_nodes, d, goal_tol, delta_min, p_polar)
    Arvore = [q_init; 0];
    convergiu = false;
    comprimento = NaN;
    
    for k = 1:max_nodes
        if rand < p_polar
            q_rand = q_goal;
        else
            q_rand = qmin + rand(7, 1) .* (qmax - qmin);
        end
        
        dists = sqrt(sum((Arvore(1:7, :) - q_rand).^2, 1));
        [~, idx_near] = min(dists);
        q_near = Arvore(1:7, idx_near);
        
        dir = q_rand - q_near;
        norm_dir = norm(dir);
        if norm_dir == 0, continue; end
        q_new = q_near + d * (dir / norm_dir);
        q_new = max(min(q_new, qmax), qmin);
        
        % Filtro de Densidade / Descarte de Nó Redundante
        dists_novas = sqrt(sum((Arvore(1:7, :) - q_new).^2, 1));
        if min(dists_novas) < delta_min
            continue;
        end
        
        if ~checa_colisao_segmento(q_near, q_new, obstaculos)
            Arvore = [Arvore, [q_new; idx_near]];
            if norm(q_new - q_goal) < goal_tol
                convergiu = true;
                break;
            end
        end
    end
    
    num_nos = size(Arvore, 2);
    if convergiu
        comprimento = calcula_comprimento_caminho(Arvore);
    end
end

%% =======================================================================
%% FUNÇÕES AUXILIARES DE COLISÃO E CINEMÁTICA
%% =======================================================================

function [q_init, q_goal] = gerar_par_valido(qmin, qmax, obstaculos)
    while true
        q_init = qmin + rand(7, 1) .* (qmax - qmin);
        if ~checa_colisao_config(q_init, obstaculos)
            break;
        end
    end
    while true
        q_goal = qmin + rand(7, 1) .* (qmax - qmin);
        if norm(q_goal - q_init) > 0.8 && ~checa_colisao_config(q_goal, obstaculos)
            break;
        end
    end
end

function colidiu = checa_colisao_segmento(q1, q2, obstaculos)
    M = 8; % Subdivisões de checagem ao longo do elo
    colidiu = false;
    for m = 0:M
        s = m / M;
        q_interp = (1 - s)*q1 + s*q2;
        if checa_colisao_config(q_interp, obstaculos)
            colidiu = true;
            return;
        end
    end
end

function colidiu = checa_colisao_config(q, obstaculos)
    colidiu = false;
    T7 = fkin_pioneer(q);
    p_efetuador = T7(1:3, 4);
    
    for i = 1:length(obstaculos)
        obs = obstaculos{i};
        if strcmp(obs.tipo, 'sphere')
            if norm(p_efetuador - obs.center) <= obs.radius
                colidiu = true; return;
            end
        elseif strcmp(obs.tipo, 'box')
            b = obs.bounds;
            if (p_efetuador(1) >= b(1,1) && p_efetuador(1) <= b(1,2)) && ...
               (p_efetuador(2) >= b(2,1) && p_efetuador(2) <= b(2,2)) && ...
               (p_efetuador(3) >= b(3,1) && p_efetuador(3) <= b(3,2))
                colidiu = true; return;
            end
        end
    end
end

function comp = calcula_comprimento_caminho(Arvore)
    caminho = Arvore(:, end);
    comp = 0;
    while true
        pai = caminho(8, end);
        if pai == 0, break; end
        q_atual = caminho(1:7, end);
        q_pai   = Arvore(1:7, pai);
        comp = comp + norm(q_atual - q_pai);
        caminho = [caminho, Arvore(:, pai)];
    end
end

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

function imprimir_metricas_pp(nome, dados)
    n_conv = sum(dados.convergidos);
    tc = (n_conv / length(dados.convergidos)) * 100;
    
    if n_conv > 0
        t_medio_ms = mean(dados.tempos(dados.convergidos)) * 1000;
        nos_m = mean(dados.nos_totais(dados.convergidos));
        comp_m = mean(dados.comprimento(dados.convergidos));
    else
        t_medio_ms = NaN; nos_m = NaN; comp_m = NaN;
    end
    
    fprintf('%-18s -> tc: %6.2f%% | t_médio: %7.2f ms | Nós na Árvore: %6.1f | Custo Caminho: %5.2f rad\n', ...
            nome, tc, t_medio_ms, nos_m, comp_m);
end