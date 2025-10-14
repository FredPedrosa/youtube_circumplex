library(readxl)
library(EGAnet)
library(readr)
library(dplyr)

data <- read_csv("~/Downloads/embeddings_circumplex.csv")

cat("Dimensões do arquivo original (Linhas x Colunas):", dim(data), "\n")
print("Primeiras linhas do arquivo original:")
print(head(data))

# --- Passo 3: Preparar o DataFrame para Análise ---
# Seleciona todas as colunas, EXCETO a primeira coluna 'palavra'
# A função select() do dplyr é mais explícita e segura para isso
embeddings_matrix <- data %>%
  select(-palavra)

print("\nÚltimas linhas da matriz de embeddings pronta para análise:")
print(tail(embeddings_matrix))

# --- Passo 1: Transpor a Matriz de Embeddings ---
transposed_matrix <- t(embeddings_matrix)

# --- Passo 2: Atribuir os Nomes das Palavras como Nomes de Coluna ---
colnames(transposed_matrix) <- data$palavra


pca_results <- prcomp(transposed_matrix, scale. = TRUE)

# --- Passo 7: Criando e Analisando o Scree Plot ---
# O scree plot nos ajuda a decidir quantas dimensões reter
screeplot(pca_results, type = "lines", npcs = 15, main = "Scree Plot dos Embeddings")
abline(h = 1, col = "red", lty = 2)

EFA.dimensions::DIMTESTS(transposed_matrix, corkind = "kendall")	

library(psych)
parallel_analysis_results <- fa.parallel(
  transposed_matrix,
  fa = "both",  # Vamos focar na Análise de Componentes Principais primeiro
  n.iter = 100, # Número de simulações
  show.legend = TRUE,
  main = "Análise Paralela de Horn (PCA)"
)



# EGA
ega_result <- EGA(data = transposed_matrix, plot.EGA = TRUE)
ega_result_TMFG <- EGA(data = transposed_matrix, model = "TMFG", plot.EGA = TRUE)

hiera <- hierEGA(data = transposed_matrix,
                 scores = "network",
                 plot.EGA = TRUE)
plot(hiera, plot.type = "separate")

#Iteração 1 - remoção de itens redundantes (UVA) e instáveis (bootEGA)
data_iter1 <- transposed_matrix
uva_results_iter1 <- UVA(data = data_iter1, verbose = FALSE)
data_after_uva_iter1 <- uva_results_iter1$reduced_data
bootEGA_results_iter1 <- bootEGA(data = data_after_uva_iter1, iter = 500, type = "parametric", model = "TMFG", seed = 123)
item_stability_proportions_iter1 <- bootEGA_results_iter1$stability$item.stability$item.stability$empirical.dimensions
stable_items_iter1 <- names(item_stability_proportions_iter1[item_stability_proportions_iter1 >= stability_threshold])
data_iter2 <- data_iter1[, stable_items_iter1]

# Iteração 2 - remoção de itens  instáveis (bootEGA)
uva_results_iter2 <- UVA(data = data_iter2, verbose = FALSE)
data_after_uva_iter2 <- uva_results_iter2$reduced_data
bootEGA_results_iter2 <- bootEGA(data = data_after_uva_iter2, iter = 500, type = "parametric", model = "TMFG", seed = 123)
item_stability_proportions_iter2 <- bootEGA_results_iter2$stability$item.stability$item.stability$empirical.dimensions
stable_items_iter2 <- names(item_stability_proportions_iter2[item_stability_proportions_iter2 >= stability_threshold])
data_iter3 <- data_iter2[, stable_items_iter2]

uva_results_iter3 <- UVA(data = data_iter3, verbose = FALSE)
data_after_uva_iter3 <- uva_results_iter3$reduced_data
bootEGA_results_iter3 <- bootEGA(data = data_after_uva_iter3, iter = 500, type = "parametric", model = "glasso", seed = 123)
item_stability_proportions_iter3 <- bootEGA_results_iter3$stability$item.stability$item.stability$empirical.dimensions
stable_items_iter3 <- names(item_stability_proportions_iter3[item_stability_proportions_iter3 >= stability_threshold])
data_iter4 <- data_iter3[, stable_items_iter3]



# verificando redundância
uva_results_iter3 <- UVA(data = data_iter3, verbose = FALSE) 



# --- Pré-requisito: Pacote lavaan e a matriz 'transposed_matrix' ---
library(lavaan)

# --- Definindo os Fatores com Base nas Palavras ---
# Fator1: Valência (Palavras positivas vs. negativas)
# Fator2: Ativação (Palavras agitadas vs. calmas)

names(transposed_matrix) <- gsub("dançando", "dancando", names(transposed_matrix))





#--- Pré-requisitos ---
library(lavaan)
# Certifique-se de que 'transposed_matrix' está carregada e 'dancando' foi renomeado.

# --- Modelo de 4 Fatores (Polos) ---
# Cada fator é definido por um conjunto de palavras o mais "puro" possível.
# Removi palavras ambíguas como 'amor' que podem ter alta ativação e alta valência.
# O objetivo aqui é fazer o modelo CONVERGIR primeiro.

modelo_4_polos_depurado <- '
  # Polo de Alta Valência e Baixa Ativação (Serenidade)
  Val_Pos_Ati_Baixa =~ paz + calma + relaxar + suave + descanse + tranquilidade
  
  # Polo de Alta Valência e Alta Ativação (Euforia)
  Val_Pos_Ati_Alta  =~ linda + feliz + amei + maravilha + perfeita + delicia + top
  
  # Polo de Baixa Valência e Alta Ativação (Tensão/Raiva)
  Val_Neg_Ati_Alta  =~ foda + merda + ruim + lixo + inferno + pedrada + pancada
  
  # Polo de Baixa Valência e Baixa Ativação (Melancolia)
  Val_Neg_Ati_Baixa =~ triste + tristeza + sozinho + chorando + chore
'

# --- Rodando a AFC ---
# Usamos `check.gradient = FALSE` e `baseline = FALSE` para ajudar na convergência inicial
# Se ainda não convergir, podemos precisar simplificar mais.
fit_4_polos <- cfa(
  modelo_4_polos_depurado, 
  data = transposed_matrix, 
  estimator = "MLR",
  verbose = TRUE, # Mostra mais detalhes do processo de otimização
  check.gradient = TRUE # Relaxa uma das checagens de convergência
)


# --- Verificando os Resultados ---
if(lavInspect(fit_4_polos, "converged")) {
  print("O MODELO CONVERGIU! Analisando os resultados...")
  
  # Pedir os índices de ajuste e os parâmetros
  summary(fit_4_polos, fit.measures = TRUE, standardized = TRUE)
  
  # Olhe especialmente a matriz de correlação entre os fatores latentes.
  # Ela nos dirá se a estrutura do circumplexo faz sentido.
  # Ex: Esperamos correlação negativa entre Val_Pos e Val_Neg.
  # E correlação próxima de zero entre Val_Pos_Ati_Baixa e Val_Neg_Ati_Alta (polos opostos).
  
} else {
  print("O MODELO NÃO CONVERGIU. Tente simplificar ainda mais os fatores.")
  print("Verifique o output do 'verbose=TRUE' para pistas.")
}


fitMeasures(fit_4_polos, fit.measures = c("chisq","df","cfi", "rmsea",
                                          "rmsea.ci.lower", "rmsea.ci.upper"))

modelo_bifatorial_final <- '
  # Fatores de grupo (os mesmos de antes)
  Serenidade =~ paz + calma + relaxar + suave + descanse + 0*tranquilidade
  Euforia   =~ linda + 0*feliz + amei + maravilha + perfeita + delicia + top
  Tensao    =~ foda + merda + ruim + lixo + inferno + pedrada + pancada
  Melancolia =~ triste + tristeza + sozinho + chorando + chore
  
  # Fator Geral, influenciando TODAS as palavras
  Geral =~ paz + calma + relaxar + suave + descanse + tranquilidade + 
           linda + feliz + amei + maravilha + perfeita + delicia + top +
           foda + merda + ruim + lixo + inferno + pedrada + pancada +
           triste + tristeza + sozinho + chorando + chore
'

# Rodar a AFC com o modelo bifatorial e ortogonalidade
fit_bifatorial <- cfa(
  modelo_bifatorial_final, 
  data = transposed_matrix, 
  estimator = "ML",
  orthogonal = TRUE # Fator Geral e fatores de grupo são independentes
)

fitMeasures(fit_bifatorial, fit.measures = c("chisq","df","cfi", "rmsea", 
                                                      "rmsea.ci.lower", "rmsea.ci.upper"))
                                             
summary(fit_bifatorial, fit.measures = TRUE, standardized = TRUE)
semTools::compRelSEM(fit_bifatorial)








# --- Pré-requisitos ---
# Pacotes 'psych' e 'GPArotation' instalados e carregados
# Matriz 'transposed_matrix' (768x45) pronta
# install.packages(c("psych", "GPArotation"))
library(psych)
library(GPArotation)

# --- Passo 1: Determinar o Número de Fatores com Análise Paralela (para EFA) ---
cat("Executando Análise Paralela para Análise Fatorial (EFA)...\n")

parallel_analysis_fa <- fa.parallel(
  transposed_matrix,
  fa = "fa", # Especifica que é para Análise Fatorial
  fm = "minres", # Método de extração (Minimum Residual - um bom padrão)
  n.iter = 100,
  main = "Análise Paralela de Horn (EFA)"
)

# A análise irá sugerir o número de fatores. Vamos supor que ela sugira 'n_fatores'.
# Pelo seu output anterior, a sugestão foi 6, mas vamos confirmar.
n_fatores_sugerido <- parallel_analysis_fa$nfact
cat("\nAnálise Paralela sugere:", n_fatores_sugerido, "fatores.\n")


# --- Passo 2: Rodar a Análise Fatorial Exploratória (EFA) ---
cat("\nRodando a EFA com", n_fatores_sugerido, "fatores e rotação oblimin...\n")

# A função fa() do pacote psych é excelente para isso.
# fm = "minres": Método de extração.
# rotate = "oblimin": Rotação oblíqua, permite que os fatores sejam correlacionados.
# nfactors: O número de fatores que queremos extrair.
efa_results <- fa(
  r = transposed_matrix,
  nfactors = n_fatores_sugerido,
  rotate = "oblimin",
  fm = "minres"
)


# --- Passo 3: Interpretar os Resultados da EFA ---

# A parte mais importante é a matriz de cargas fatoriais (Factor Loadings).
# A função print() do psych já formata isso de uma maneira muito útil.
# Ela esconde cargas baixas (por padrão, < 0.3) para facilitar a visualização.
print(efa_results$loadings, cutoff = 0.3, sort = TRUE)

# Analise a tabela de cargas:
# - Cada coluna (MR1, MR2, etc.) é um fator.
# - Olhe quais palavras têm cargas altas em cada fator.
# - Tente "nomear" cada fator com base no tema comum das palavras que carregam nele.
#   - Ex: Se "triste", "sozinho", "chore" carregam alto no Fator 1, você pode chamá-lo de "Melancolia".
#   - Se "paz", "calma", "relaxar" carregam alto no Fator 2, é "Serenidade".

# Veja também a correlação entre os fatores
cat("\nMatriz de Correlação entre os Fatores (Factor Correlation):\n")
print(round(efa_results$Phi, 2))
# Isso lhe dirá como os seus fatores recém-descobertos se relacionam.





# --- PASSO 1: PURIFICAÇÃO ---
palavras_a_remover <- c("pancada", "inferno", "recordações", "bons", "sozinho", 
                        "suave", "amor", "lembro", "top", "viciada", "chore", "graça", "dispara")

# Mantém apenas as colunas (palavras) que NÃO estão na lista de remoção
matriz_purificada <- transposed_matrix[, !(colnames(transposed_matrix) %in% palavras_a_remover)]

# --- PASSO 2: NOVA EFA ---
parallel_analysis_fa <- fa.parallel(
  matriz_purificada,
  fa = "fa", # Especifica que é para Análise Fatorial
  fm = "minres", # Método de extração (Minimum Residual - um bom padrão)
  n.iter = 100,
  main = "Análise Paralela de Horn (EFA)"
)

efa_results <- fa(
  r = transposed_matrix,
  nfactors = 4,
  rotate = "oblimin",
  fm = "minres"
  )
print(efa_results$loadings, cutoff = 0.3, sort = TRUE)


  


# --- Passo 3: Interpretar os Resultados da EFA ---

# A parte mais importante é a matriz de cargas fatoriais (Factor Loadings).
# A função print() do psych já formata isso de uma maneira muito útil.
# Ela esconde cargas baixas (por padrão, < 0.3) para facilitar a visualização.
print(efa_results$loadings, cutoff = 0.3, sort = TRUE)

# Analise a tabela de cargas:
# - Cada coluna (MR1, MR2, etc.) é um fator.
# - Olhe quais palavras têm cargas altas em cada fator.
# - Tente "nomear" cada fator com base no tema comum das palavras que carregam nele.
#   - Ex: Se "triste", "sozinho", "chore" carregam alto no Fator 1, você pode chamá-lo de "Melancolia".
#   - Se "paz", "calma", "relaxar" carregam alto no Fator 2, é "Serenidade".

# Veja também a correlação entre os fatores
cat("\nMatriz de Correlação entre os Fatores (Factor Correlation):\n")
print(round(efa_results$Phi, 2))







palavras_a_remover_2 <- c(
  "grande", "bons", "vibe", "suave", "amei", "gosto", "viciada",
  "inferno", "top", "delicia", 
  "tranquilidade", # Idem
  "paz",    # Idem
  "feliz", # a partir daqui foi a 3a rodada 
  "relaxar",
  "recordações",
  "pancada",
  "maravilha",
  "amor", # 4a rodada
  "lembro"
)

# Mantém apenas as colunas (palavras) que NÃO estão na lista de remoção
matriz_purificada_2 <- transposed_matrix[, !(colnames(transposed_matrix) %in% palavras_a_remover_2)]

cat("Número de palavras restantes para a análise:", ncol(matriz_purificada_2), "\n")


# --- PASSO 2: NOVA ANÁLISE PARALELA ---
parallel_analysis_fa_2 <- fa.parallel(
  matriz_purificada_2,
  fa = "fa",
  fm = "minres",
  n.iter = 100,
  main = "Análise Paralela (Rodada 2)"
)

n_fatores_sugerido_2 <- parallel_analysis_fa_2$nfact
cat("\nAnálise Paralela na matriz purificada sugere:", n_fatores_sugerido_2, "fatores.\n")


# --- PASSO 3: NOVA EFA - TESTANDO 2 E 3 FATORES ---

# Teste com 4 Fatores 
cat("\n--- EFA com 4 Fatores ---\n")
efa_4_fatores <- fa(
  r = matriz_purificada_2,
  nfactors = 4,
  rotate = "oblimin",
  fm = "minres"
)
print(efa_4_fatores$loadings, cutoff = 0.3, sort = TRUE)
cat("\nCorrelações entre 2 Fatores:\n")
print(round(efa_2_fatores$Phi, 2))



# --- Pré-requisito: Matriz 'transposed_matrix' (768x45) e a lista 'palavras_a_remover_2' ---

# --- PASSO DE PURIFICAÇÃO (Já feito por você) ---
matriz_purificada_2 <- transposed_matrix[, !(colnames(transposed_matrix) %in% palavras_a_remover_2)]


# --- NOVO PASSO: SANITIZAR NOMES DE COLUNAS ---
# Esta função irá remover todos os acentos e caracteres especiais dos nomes das colunas

# install.packages("stringi") # Se não tiver instalado
library(stringi)

# Remove acentos e caracteres especiais dos nomes de todas as colunas
colnames(matriz_purificada_2) <- stri_trans_general(colnames(matriz_purificada_2), "Latin-ASCII")

# Verifica os novos nomes para garantir que 'graca' e 'coracao' foram corrigidos
print("Nomes de colunas sanitizados:")
print(colnames(matriz_purificada_2))


# --- AGORA, A SINTAXE DA AFC IRÁ FUNCIONAR ---
library(lavaan)

# --- MODELO 1: AFC com 4 Fatores Correlacionados (Sintaxe Corrigida) ---
cat("\n--- Rodando Modelo 1: 4 Fatores Correlacionados ---\n")

modelo_4_fatores_final <- '
  # Fatores de Primeira Ordem (Polos do Circumplexo)
  Quietude    =~ saudade + descanse + triste + sozinho + tristeza + chore + chorando + coracao + graca + dispara + penso + calma
  Apreciacao  =~ linda + maravilhosa + amo + perfeita + gostei
  Aversao     =~ foda + merda + ruim + lixo + pedrada
  Exuberancia =~ louco + dancando + charmosa + chique
'

fit_4_fatores <- cfa(modelo_4_fatores_final, data = matriz_purificada_2, estimator = "MLR")
print("--- Resultados Modelo 1 (4 Fatores) ---")
summary(fit_4_fatores, fit.measures = TRUE, standardized = TRUE)


# --- MODELO 2: AFC de Segunda Ordem (Sintaxe Corrigida) ---
cat("\n--- Rodando Modelo 2: Fator de Segunda Ordem ---\n")

modelo_2a_ordem_final <- '
  # Fatores de Primeira Ordem
  Quietude    =~ saudade + descanse + triste + sozinho + tristeza + chore + chorando + coracao + graca + dispara + penso + calma
  Apreciacao  =~ linda + maravilhosa + amo + perfeita + gostei
  Aversao     =~ foda + merda + ruim + lixo + pedrada
  Exuberancia =~ louco + dancando + charmosa + chique
  
  # Fator de Segunda Ordem
  EmocionalidadeGeral =~ Quietude + Apreciacao + Aversao + Exuberancia
'

fit_2a_ordem <- cfa(modelo_2a_ordem_final, data = matriz_purificada_2, estimator = "MLR")
print("--- Resultados Modelo 2 (Segunda Ordem) ---")
summary(fit_2a_ordem, fit.measures = TRUE, standardized = TRUE)


# --- MODELO 3: AFC Bifatorial (Sintaxe Corrigida) ---
cat("\n--- Rodando Modelo 3: Modelo Bifatorial ---\n")

# A função 'make.names' é outra forma de garantir nomes válidos
todas_palavras_purificadas_sanitizadas <- paste(make.names(colnames(matriz_purificada_2)), collapse = " + ")

modelo_bifatorial_final <- '
  # Fatores de Grupo Específicos
  Quietude    =~ saudade + descanse + triste + sozinho + tristeza + chore + chorando + coracao + graca + dispara + penso + calma
  Apreciacao  =~ linda + maravilhosa + amo + perfeita + gostei
  Aversao     =~ foda + merda + ruim + lixo + pedrada
  Exuberancia =~ louco + dancando + charmosa + chique

  # Fator Geral
  Geral =~ saudade + descanse + triste + sozinho + tristeza + chore + chorando + coracao + graca + dispara + penso + calma +
           linda + maravilhosa + amo + perfeita + gostei +
           foda + merda + ruim + lixo + pedrada +
           louco + dancando + charmosa + chique
'

fit_bifatorial <- cfa(modelo_bifatorial_final, 
                      data = matriz_purificada_2, 
                      estimator = "MLR", 
                      orthogonal = TRUE)
print("--- Resultados Modelo 3 (Bifatorial) ---")
summary(fit_bifatorial, fit.measures = TRUE, standardized = TRUE)


# --- PASSO FINAL: COMPARAR OS MODELOS ---
cat("\n--- Comparando os Índices de Ajuste dos Modelos ---\n")
anova(fit_4_fatores, fit_2a_ordem, fit_bifatorial)

























# --- Pré-requisito: Pacote 'lavaan' e a matriz 'matriz_purificada_2' ---
library(lavaan)

# --- MODELO 1: AFC com 4 Fatores Correlacionados (Sintaxe Correta) ---
cat("\n--- Rodando Modelo 1: 4 Fatores Correlacionados ---\n")

modelo_4_fatores_final <- '
  # Fatores de Primeira Ordem (Polos do Circumplexo)
  Quietude    =~ saudade + descanse + triste + sozinho + tristeza + chore + chorando + coração + graça + dispara + penso + calma
  Apreciacao  =~ linda + maravilhosa + amo + perfeita + gostei
  Aversao     =~ foda + merda + ruim + lixo + pedrada
  Exuberancia =~ louco + dancando + charmosa + chique
'

fit_4_fatores <- cfa(modelo_4_fatores_final, data = matriz_purificada_2, estimator = "MLR")
print("--- Resultados Modelo 1 (4 Fatores) ---")
summary(fit_4_fatores, fit.measures = TRUE, standardized = TRUE)


# --- MODELO 2: AFC de Segunda Ordem (Sintaxe Correta) ---
cat("\n--- Rodando Modelo 2: Fator de Segunda Ordem ---\n")

modelo_2a_ordem_final <- '
  # Fatores de Primeira Ordem
  Quietude    =~ saudade + descanse + triste + sozinho + tristeza + chore + chorando + coração + graça + dispara + penso + calma
  Apreciacao  =~ linda + maravilhosa + amo + perfeita + gostei
  Aversao     =~ foda + merda + ruim + lixo + pedrada
  Exuberancia =~ louco + dancando + charmosa + chique
  
  # Fator de Segunda Ordem
  EmocionalidadeGeral =~ Quietude + Apreciacao + Aversao + Exuberancia
'

fit_2a_ordem <- cfa(modelo_2a_ordem_final, data = matriz_purificada_2, estimator = "MLR")
print("--- Resultados Modelo 2 (Segunda Ordem) ---")
summary(fit_2a_ordem, fit.measures = TRUE, standardized = TRUE)


# --- MODELO 3: AFC Bifatorial (Sintaxe Correta) ---
cat("\n--- Rodando Modelo 3: Modelo Bifatorial ---\n")

modelo_bifatorial_final <- '
  # Fatores de Grupo Específicos
  Quietude    =~ saudade + descanse + triste + sozinho + tristeza + chore + chorando + coracao + graca + dispara + penso + calma
  Apreciacao  =~ linda + maravilhosa + amo + perfeita + gostei
  Aversao     =~ foda + merda + ruim + lixo + pedrada
  Exuberancia =~ louco + dancando + charmosa + chique

  # Fator Geral
  Geral =~ saudade + descanse + triste + sozinho + tristeza + chore + chorando + coracao + graca + dispara + penso + calma +
           linda + maravilhosa + amo + perfeita + gostei +
           foda + merda + ruim + lixo + pedrada +
           louco + dancando + charmosa + chique
'

fit_bifatorial <- cfa(modelo_bifatorial_final, 
                      data = matriz_purificada_2, 
                      estimator = "MLR", 
                      orthogonal = TRUE) # Crucial para o modelo bifatorial
print("--- Resultados Modelo 3 (Bifatorial) ---")
summary(fit_bifatorial, fit.measures = TRUE, standardized = TRUE)


# --- PASSO FINAL: COMPARAR OS MODELOS ---
cat("\n--- Comparando os Índices de Ajuste dos Modelos ---\n")
# A função anova() compara modelos aninhados. O modelo de 4 fatores é aninhado nos outros.
# O bifatorial e o de segunda ordem não são aninhados entre si, então comparamos seus AIC/BIC.
anova(fit_4_fatores, fit_2a_ordem, fit_bifatorial)

# Para comparar os modelos não aninhados (2a ordem vs. Bifatorial), olhamos o AIC e BIC
# no output de cada summary(). O menor valor é preferível.













# --- Pré-requisitos ---
# Pacotes 'psych' e 'ggplot2' (para gráficos melhores)
# Matriz 'transposed_matrix' (768x45) pronta
# install.packages(c("psych", "ggplot2", "ggrepel"))
library(psych)
library(ggplot2)
library(ggrepel)

# --- PASSO 1: Rodar a Análise de Componentes Principais ---
# Usando a função principal() do pacote psych
# nfactors = 3: Extraímos 3 componentes, conforme sugerido pela Análise Paralela
# rotate = "none": É crucial NÃO rotacionar inicialmente, para ver o fator geral.

cat("\n--- Rodando PCA para extrair 3 componentes ---\n")
pca_results_psych <- principal(
  r = transposed_matrix,
  nfactors = 3,
  rotate = "none" # Sem rotação para ver a estrutura bruta
)

# --- PASSO 2: Analisar a Matriz de Cargas ---
print(pca_results_psych$loadings, cutoff = 0.3, sort = TRUE)
pca_results_psych$fit.off

# Interpretação esperada:
# - PC1: Quase todas as palavras devem ter cargas altas e positivas (Fator Geral).
# - PC2 e PC3: As cargas serão mistas (positivas e negativas) e devem conter a estrutura
#   dos polos do circumplexo.

# --- PASSO 3: Criar o Gráfico de Dispersão dos Componentes 2 e 3 ---
# Primeiro, extraímos as cargas em um dataframe para usar com ggplot2
loadings_df <- as.data.frame(unclass(pca_results_psych$loadings))

# Adicionar os nomes das palavras como uma coluna
loadings_df$palavra <- rownames(loadings_df)

cat("\n--- Gerando o mapa afetivo com os Componentes 2 e 3 ---\n")

# Criar o gráfico com ggplot2
pca_plot <- ggplot(loadings_df, aes(x = PC2, y = PC3, label = palavra)) +
  # Adiciona os pontos
  geom_point(color = "red", size = 2) +
  
  # Adiciona os eixos x=0 e y=0 para formar os quadrantes
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # Adiciona os rótulos das palavras de forma inteligente para evitar sobreposição
  geom_text_repel(size = 3.5, max.overlaps = Inf) +
  
  # Títulos e labels
  labs(
    title = "Mapa Afetivo a partir dos Componentes Principais 2 e 3",
    subtitle = "Estrutura após a extração do Fator Geral (PC1)",
    x = paste0("Componente Principal 2 (", round(pca_results_psych$Vaccounted[2, "Proportion Var"] * 100, 1), "% da Variância)"),
    y = paste0("Componente Principal 3 (", round(pca_results_psych$Vaccounted[3, "Proportion Var"] * 100, 1), "% da Variância)")
  ) +
  
  # Tema limpo
  theme_minimal() +
  
  # Garante que os eixos tenham a mesma escala para uma visualização correta do círculo
  coord_fixed()

# Exibir o gráfico
print(pca_plot)




# --- PASSO 3: Criar o Gráfico de Dispersão dos Componentes 2 e 3 (CORRIGIDO) ---
# Primeiro, extraímos as cargas em um dataframe para usar com ggplot2
loadings_df <- as.data.frame(unclass(pca_results_psych$loadings))

# Adicionar os nomes das palavras como uma coluna
loadings_df$palavra <- rownames(loadings_df)

cat("\n--- Gerando o mapa afetivo com os Componentes 2 e 3 ---\n")

# Extrair a variância explicada de forma correta
variancia_explicada_pca <- pca_results_psych$Vaccounted["Proportion Var",]

# Criar o gráfico com ggplot2
pca_plot <- ggplot(loadings_df, aes(x = PC2, y = PC3, label = palavra)) +
  geom_point(color = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text_repel(size = 3.5, max.overlaps = Inf) +
  labs(
    title = "Mapa Afetivo a partir dos Componentes Principais 2 e 3",
    subtitle = "Estrutura após a extração do Fator Geral (PC1)",
    x = paste0("Componente Principal 2 (", round(variancia_explicada_pca['PC2'] * 100, 1), "% da Variância)"),
    y = paste0("Componente Principal 3 (", round(variancia_explicada_pca['PC3'] * 100, 1), "% da Variância)")
  ) +
  theme_minimal() +
  coord_fixed()

# Exibir o gráfico
print(pca_plot)




# --- Pré-requisitos ---
# Objeto 'pca_results_psych' e 'loadings_df' já criados do código anterior.
# Pacotes 'ggplot2' e 'ggrepel' carregados.

# --- Gerando o Gráfico 2D com a Terceira Dimensão (Intensidade) ---

cat("\n--- Gerando o mapa afetivo 2D com a dimensão de Intensidade (PC1) ---\n")

# O dataframe 'loadings_df' já contém as colunas PC1, PC2, PC3 e palavra.

plot_2d_com_intensidade <- ggplot(loadings_df, aes(x = PC3, y = PC2, label = palavra)) +
  # Adiciona os eixos x=0 e y=0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray10") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray10") +
  
  # AQUI ESTÁ A MÁGICA:
  # Mapeia a carga no PC1 para a COR e o TAMANHO dos pontos.
  # Pontos com alta carga em PC1 (alta intensidade) serão maiores e mais escuros.
  geom_point(aes(color = PC1, size = PC1), alpha = 0.7) +
  
  # Adiciona os rótulos das palavras
  geom_text_repel(size = 2.5, max.overlaps = Inf, segment.color = "grey50") +
  
  # Define a paleta de cores (ex: de azul claro para azul escuro)
  scale_color_gradient(low = "purple", high = "orange") +
  
  # Define a escala de tamanho
  scale_size(range = c(2, 10)) + # Pontos variarão de tamanho 2 a 10
  
  # Títulos e labels
  labs(
    title = "Mapa Afetivo com Dimensão de Intensidade",
    subtitle = "Circumplexo (PC2 vs PC3) com Intensidade (PC1) codificada na cor e tamanho",
    x = "Componente 3 (Valência Invertida)",
    y = "Componente 2 (Ativação)",
    color = "Intensidade (Carga PC1)", # Título da legenda de cor
    size = "Intensidade (Carga PC1)"  # Título da legenda de tamanho
  ) +
  
  theme_minimal() +
  coord_fixed()

# Exibir o gráfico
print(plot_2d_com_intensidade)





# --- Pré-requisitos ---
# Pacotes 'dplyr', 'ggplot2', 'ggrepel', 'readr'
# Objeto 'loadings_df' já criado com as colunas PC1, PC2, PC3, e palavra.

library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)

# --- PASSO 1: Adicionar a Frequência das Palavras ao Dataframe de Plotagem ---

# Carregar o arquivo que contém as frequências
# (Mesmo que você já tenha carregado antes, vamos fazer de novo para garantir)
# Lembre-se de ajustar o caminho se necessário
freq_data <- read_csv("~/Downloads/tabela_top10_palavras_por_fator.csv")

# Renomear as colunas para facilitar a junção
# O nome da coluna pode variar dependendo de como você o salvou. Verifique o seu CSV!
# Vamos supor que os nomes sejam 'Palavra-Chave Afetiva' e 'Frequência de Ocorrência'
freq_data <- freq_data %>%
  rename(
    palavra = `Palavra-Chave Afetiva`,
    frequencia = `Frequência de Ocorrência`
  ) %>%
  select(palavra, frequencia) %>% # Seleciona apenas as colunas que nos interessam
  distinct(palavra, .keep_all = TRUE) # Garante que cada palavra apareça apenas uma vez

# Agora, vamos juntar (merge) as frequências ao nosso dataframe de cargas fatoriais
loadings_df_com_freq <- left_join(loadings_df, freq_data, by = "palavra")

# Verificar se a junção funcionou
print(head(loadings_df_com_freq))


# --- PASSO 2: Encontrar os limites para o gráfico quadrado ---
max_limit <- max(abs(c(loadings_df_com_freq$PC2, loadings_df_com_freq$PC3))) * 1.1


# --- PASSO 3: Gerar o Gráfico Final (Cor = Intensidade, Tamanho = Frequência) ---

plot_final_informativo <- ggplot(loadings_df_com_freq, aes(x = PC3, y = PC2, label = palavra)) +
  # Eixos
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  
  # AQUI ESTÁ A NOVA MÁGICA:
  # color = PC1 (Intensidade)
  # size = frequencia (Frequência de Uso)
  geom_point(aes(color = PC1, size = frequencia), alpha = 0.7) +
  
  # Rótulos
  geom_text_repel(
    size = 3.5, 
    max.overlaps = Inf, 
    segment.color = "grey50",
    point.padding = 0.5, 
    box.padding = 0.5
  ) +
  
  # Escala de cores para Intensidade
  scale_color_gradient(low = "purple", high = "orange") +
  
  # Escala de tamanho para Frequência
  scale_size(range = c(3, 12)) + # Ajuste os valores min e max para o melhor efeito visual
  
  # Limites dos eixos para ser quadrado
  xlim(-max_limit, max_limit) +
  ylim(-max_limit, max_limit) +
  
  # Títulos e legendas
  labs(
    title = "Mapa Afetivo Hierárquico",
    subtitle = "Circumplexo (Posição), Intensidade (Cor) e Frequência (Tamanho)",
    x = "Componente 3 (Valência Invertida)",
    y = "Componente 2 (Ativação)",
    color = "Intensidade\n(Carga PC1)",
    size = "Frequência\nde Uso"
  ) +
  
  # Mantém a escala 1:1
  coord_fixed() +
  
  # Tema
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

# Exibir o gráfico
print(plot_final_informativo)


# --- Pré-requisitos ---
# Objeto 'loadings_df' já criado.
# Pacote 'plotly' instalado e carregado.
# install.packages("plotly")
library(plotly)

cat("\n--- Gerando o mapa afetivo 3D interativo (PC1, PC2, PC3) ---\n")

# Mapeando os eixos para a nossa interpretação teórica
# Eixo Z: PC1 (Intensidade)
# Eixo X: PC2 (Ativação)
# Eixo Y: PC3 (Valência Invertida)

fig_3d <- plot_ly(
  data = loadings_df, 
  x = ~PC2, 
  y = ~PC3, 
  z = ~PC1, 
  type = 'scatter3d', 
  mode = 'markers+text',
  text = ~palavra, # Mostra a palavra ao lado do ponto
  textfont = list(size = 10),
  marker = list(
    size = 5,
    color = ~PC1, # Colore os pontos pela sua altura no eixo Z (Intensidade)
    colorscale = 'Viridis', # Uma boa paleta de cores
    showscale = TRUE,
    colorbar = list(title = "Intensidade (PC1)")
  )
)

# Adiciona títulos e labels aos eixos
fig_3d <- fig_3d %>% layout(
  title = "Estrutura Hierárquica do Afeto (Visualização 3D)",
  scene = list(
    xaxis = list(title = "Ativação (PC2)"),
    yaxis = list(title = "Valência Invertida (PC3)"),
    zaxis = list(title = "Intensidade Geral (PC1)")
  )
)

# Exibir o gráfico interativo
fig_3d





#### NÃO FUNCIONA

library(seminr)

# Seus dados transpostos, onde as colunas são as palavras
# transposed_matrix <- ...

# Converter para data.frame
dados_para_seminr <- as.data.frame(transposed_matrix)

# Modelo de Medição: Um único fator formativo gigante
# que inclui todas as palavras.
measurement_model_formativo <- constructs(
  composite("Intensidade", names(dados_para_seminr), weights = mode_B)
)

# Modelo Estrutural Vazio
structural_model_vazio <- relationships(
  paths(from = "Intensidade", to = c())
)

# Estimar o modelo (pode falhar por causa dos bugs que vimos)
pls_model <- estimate_pls(
  data = dados_para_seminr,
  measurement_model = measurement_model_formativo,
  structural_model = structural_model_vazio
)

# Se funcionar, o mais importante é o bootstrapping:
boot_pls <- bootstrap_model(pls_model, nboot = 500)
summary_boot <- summary(boot_pls)

# Verificar a estabilidade dos pesos
print(summary_boot$bootstrapped_weights)





#install.packages("cSEM")
library(cSEM)
library(janitor)
data <- read_csv("~/Downloads/embeddings_circumplex.csv")

# Extrair os nomes das palavras
nomes_palavras <- data$palavra

# LIMPEZA DOS NOMES DAS PALAVRAS
# A função clean_names() do janitor remove acentos, espaços e caracteres especiais
nomes_limpos <- janitor::make_clean_names(nomes_palavras)

# Diagnóstico: veja a diferença
print("--- Comparando nomes originais e limpos ---")
print(data.frame(original = nomes_palavras, limpo = nomes_limpos))

# Preparar a matriz de dados para a análise
embeddings_matrix <- data %>% select(-palavra)
transposed_matrix <- t(embeddings_matrix)

# ATRIBUIR OS NOMES LIMPOS
colnames(transposed_matrix) <- nomes_limpos

# Converter para data.frame
dados_csem <- as.data.frame(transposed_matrix)


# --------------------------------------------------------------------------
# PASSO 3: ESPECIFICAR O MODELO USANDO OS NOMES LIMPOS
# --------------------------------------------------------------------------
# O cSEM agora receberá um modelo apenas com caracteres válidos.
# O truque é construir a fórmula programaticamente para não ter que digitar tudo.

# Pega todos os nomes de colunas (já limpos)
todos_os_indicadores <- names(dados_csem)

# Constrói a string da fórmula
formula_formativa <- paste(todos_os_indicadores, collapse = " + ")

# Cria o modelo final
modelo_csem_limpo <- paste("Intensidade <~", formula_formativa)

# Diagnóstico: veja como ficou o modelo
print("\n--- Modelo CSEM com Nomes Limpos ---")
print(modelo_csem_limpo)


# --------------------------------------------------------------------------
# PASSO 4: RODAR A ANÁLISE COM cSEM
# --------------------------------------------------------------------------
print("\n--- Tentando estimar com cSEM e nomes de variáveis limpos ---")
try({
  csem_results <- csem(
    .data = dados_csem,
    .model = modelo_csem_limpo
  )
  
  print("************************************************************")
  print("FUNCIONOU! O problema eram os caracteres especiais nos nomes.")
  print("************************************************************")
  
  # Avaliar os resultados
  resumo_csem <- summary(csem_results)
  print(resumo_csem)
  
}, silent = FALSE)






# --------------------------------------------------------------------------
# PASSO 1: PREPARAR OS DADOS (JÁ FEITO, USANDO 'dados_csem' COM NOMES LIMPOS)
# --------------------------------------------------------------------------
# install.packages("cSEM")
# install.packages("janitor")
library(cSEM)
library(janitor)
library(readr)
library(dplyr)

# (Se necessário, rodar o código de limpeza de nomes de novo)
data <- read_csv("~/Downloads/embeddings_circumplex.csv")
nomes_palavras <- data$palavra
nomes_limpos <- make_clean_names(nomes_palavras)
embeddings_matrix <- data %>% select(-palavra)
transposed_matrix <- t(embeddings_matrix)
colnames(transposed_matrix) <- nomes_limpos
dados_csem <- as.data.frame(transposed_matrix)


# --------------------------------------------------------------------------
# PASSO 2: ESPECIFICAR O MODELO DE 3 FATORES FORMATIVOS
# --------------------------------------------------------------------------
# Atribuição baseada na interpretação do seu Mapa Afetivo (PCA Biplot)

modelo_3_fatores <- '
  # Measurement model (todos formativos)

  # Fator Intensidade: Todas as palavras, como sugerido pela dominância do PC1
  Intensidade <~ foda + grande + saudade + merda + bons + descanse +
                 delicia + triste + sozinho + tristeza + maravilha +
                 louco + dancando + charmosa + chique + linda + paz +
                 feliz + vibe + maravilhosa + chore + chorando + suave +
                 relaxar + tranquilidade + amo + amei + amor + gosto +
                 coracao + lembro + graca + dispara + penso + calma +
                 top + perfeita + ruim + gostei + lixo + viciada +
                 pedrada + inferno + pancada + recordacoes

  # Fator Ativação: Palavras dos polos do PC2
  Ativacao <~ dancando + charmosa + chique + louco + top + viciada + pancada +
               paz + relaxar + tranquilidade + calma + descanse + suave

  # Fator Valência: Palavras dos polos do PC3
  Valencia <~ foda + merda + ruim + lixo + inferno + triste + tristeza +
               amor + amei + amo + linda + maravilha + maravilhosa + feliz +
               gosto + gostei + bons + perfeita + paz

  # Structural model: Permitir que todos os fatores se correlacionem
  Intensidade ~~ Ativacao
  Intensidade ~~ Valencia
  Ativacao    ~~ Valencia
'

# --------------------------------------------------------------------------
# PASSO 3: RODAR A ANÁLISE COM cSEM
# --------------------------------------------------------------------------
print("\n--- Tentando estimar o modelo de 3 fatores com cSEM ---")
try({
  csem_3f_results <- csem(
    .data = dados_csem,
    .model = modelo_3_fatores
  )
  
  print("************************************************************")
  print("FUNCIONOU! O modelo de 3 fatores foi estimado.")
  print("************************************************************")
  
  # Avaliar os resultados
  resumo_csem_3f <- summary(csem_3f_results)
  print(resumo_csem_3f)
  
  # Bootstrapping para os pesos
  print("\nIniciando bootstrapping...")
  boot_csem_3f <- bootstrap(csem_3f_results, .R = 500)
  
  print("\n--- Resultados do Bootstrapping (p-valores para os pesos) ---")
  print(summary(boot_csem_3f)$Estimates$Bootstrapped_weights)
  
}, silent = FALSE)
