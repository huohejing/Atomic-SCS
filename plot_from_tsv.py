#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

def extract_scores(tsv_file):
    """从 atomic_scs.py 输出的 TSV 文件中提取分数（第二列）"""
    scores = []
    with open(tsv_file, 'r') as f:
        lines = f.readlines()
        if not lines:
            return scores
        # 第一行是标题，跳过
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    scores.append(float(parts[1]))
                except ValueError:
                    scores.append(1.0)
    return scores

# 提取分数
normal_scores = extract_scores("normal_raw.tsv")
problematic_scores = extract_scores("abnormal_raw.tsv")   # 如果你的文件已改名为 problematic_raw.tsv，这里也改

print(f"正常分子数量: {len(normal_scores)}")
print(f"问题分子数量: {len(problematic_scores)}")

# 构建 DataFrame
df_normal = pd.DataFrame({'Score': normal_scores, 'Group': 'Normal'})
df_problematic = pd.DataFrame({'Score': problematic_scores, 'Group': 'Problematic'})
data = pd.concat([df_normal, df_problematic], ignore_index=True)

# 绘制箱线图
plt.figure(figsize=(6, 8))
sns.boxplot(x='Group', y='Score', data=data, palette=['#1f77b4', '#ff7f0e'])
sns.stripplot(x='Group', y='Score', data=data, color='black', alpha=0.3, size=2)
plt.title('Atomic-SCS Scores: Normal vs. Problematic Molecules')
plt.ylabel('Compliance Score (lower is better)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('boxplot.png', dpi=300)
plt.savefig('boxplot.pdf')
print("箱线图已保存为 boxplot.png 和 boxplot.pdf")

# 统计检验（正常组分数是否显著低于异常组）
u_stat, p_value = mannwhitneyu(normal_scores, problematic_scores, alternative='less')
print(f"Mann-Whitney U 检验: U = {u_stat:.2f}, p = {p_value:.2e}")