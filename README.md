# Análise de Conectividade Funcional em EEG com Coerência (MATLAB)

Este projeto analisa a conectividade funcional em sinais EEG recolhidos em três estados distintos: tarefa de cálculo ativa com olhos fechados, repouso com olhos fechados e repouso com olhos abertos.

## Objetivo

Avaliar as diferenças de conectividade entre canais EEG através da coerência, com foco em bandas específicas (theta, alpha, beta), e representá-las com matrizes de coerência de forma a poder identificar padrões indicativos dos três estados diferentes.

## Metodologia

- Pré-processamento com remoção de offset e filtros (notch, passa-alto e passa-baixo)
- Filtragem por banda (theta: 4–8 Hz, alpha: 8–13 Hz, etc.)
- Cálculo da coerência entre todos os pares de canais com `mscohere`
- Geração de heatmaps para visualizar padrões espaciais de conectividade

## Tecnologias usadas

- MATLAB
- EEGLAB (para `topoplot` e leitura de locais dos eletrodos)
- Sinais EEG reais (.bdf)

## Conteúdo

- `main.m`: script principal que executa toda a pipeline
- `apply_filters.m`: função para aplicar múltiplos filtros a um sinal
- `calculate_coherence.m`: função que calcula a matriz de coerência entre canais
- `projeto_final_v3.pdf`: relatório final do projeto

## Exemplos de resultados

Inclui:
- Matrizes de coerência por condição
- Mapas topográficos das médias por canal
- Diferença de coerência entre estados (ex: ativo vs. repouso)

## Autor

João Mota — Estudante de Engenharia Biomédica na NOVA FCT

## Repositório

Este projeto foi desenvolvido no contexto da unidade curricular de **Eletrofisiologia**.

