import pandas as pd
from scipy.stats import spearmanr

# Load the dataset
file_name = 'EngChain_Model_Comparison.csv'
df = pd.read_csv(file_name)

# Calculate the Spearman correlation coefficient and the p-value
correlation, p_value = spearmanr(df['final_answer_match'], df['step_f1'])

print(f"Spearman Correlation Coefficient: {correlation:.4f}")
print(f"P-value: {p_value:.4e}")

# Interpretation
if p_value < 0.05:
    print("The correlation is statistically significant (p < 0.05).")
else:
    print("The correlation is not statistically significant (p >= 0.05).")
