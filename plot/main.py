import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Carregando os dados
df_residuos = pd.read_csv("../data/benchmark_residuos.csv")
df_tempos = pd.read_csv("../data/benchmark_tempos.csv")

# Configurando o estilo dos gráficos
sns.set_theme(style="whitegrid")

# 1. Gráfico de Tempo vs Tamanho
plt.figure(figsize=(10, 6))
sns.lineplot(data=df_tempos, x='Tamanho_N', y='Tempo_s', hue='Metodo', marker='o', linewidth=2)
plt.yscale('log')
plt.xscale('log')
plt.title('Tempo de Execução vs Tamanho do Sistema (N)', fontsize=14)
plt.xlabel('Tamanho N (escala log)', fontsize=12)
plt.ylabel('Tempo em segundos (escala log)', fontsize=12)
plt.legend(title='Método')
plt.tight_layout()
plt.savefig('tempo_vs_tamanho.png', dpi=300)

# 2. Gráfico de Norma do Resíduo vs Tamanho
plt.figure(figsize=(10, 6))
sns.lineplot(data=df_residuos, x='Tamanho_N', y='Norma_Residuo', hue='Metodo', marker='s', linewidth=2)
plt.yscale('log')
plt.title('Norma do Resíduo vs Tamanho do Sistema (N)', fontsize=14)
plt.xlabel('Tamanho N', fontsize=12)
plt.ylabel('Norma do Resíduo (escala log)', fontsize=12)
plt.legend(title='Método')
plt.tight_layout()
plt.savefig('residuo_vs_tamanho.png', dpi=300)

# 3. Gráfico de Trade-off (Tempo vs Precisão) para o maior N
max_n = df_tempos['Tamanho_N'].max()
df_merged = pd.merge(df_tempos[df_tempos['Tamanho_N'] == max_n],
                     df_residuos[df_residuos['Tamanho_N'] == max_n],
                     on=['Tamanho_N', 'Metodo'])

plt.figure(figsize=(10, 6))
sns.scatterplot(data=df_merged, x='Tempo_s', y='Norma_Residuo', hue='Metodo', s=300, edgecolor='black')
plt.yscale('log')
plt.xscale('log')
plt.title(f'Trade-off: Tempo vs Resíduo (para N = {max_n})', fontsize=14)
plt.xlabel('Tempo (s)', fontsize=12)
plt.ylabel('Norma do Resíduo', fontsize=12)

# Adicionando os nomes dos métodos direto no gráfico para melhor visualização
for i in range(df_merged.shape[0]):
    plt.text(df_merged['Tempo_s'].iloc[i] * 1.1, 
             df_merged['Norma_Residuo'].iloc[i], 
             df_merged['Metodo'].iloc[i], 
             fontsize=10)

plt.legend([],[], frameon=False)
plt.tight_layout()
plt.savefig('tradeoff.png', dpi=300)