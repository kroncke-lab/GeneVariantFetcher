import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Load Data
file_path = "/Users/kronckbm/GitRepos/GeneVariantFetcher/tests/TTR_Penetrance_Data.csv"
df = pd.read_csv(file_path)

# 2. Filter valid data
# Ensure we strictly look at the calculated penetrance column
data_to_plot = df["Calc_Penetrance_Pct"].dropna()

# 3. Plot Histogram
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))

# Histogram with Kernel Density Estimate (KDE) line
sns.histplot(
    data=data_to_plot,
    bins=20,  # Breaks data into 5% intervals (100/20)
    kde=True,  # Adds the smooth curve
    color="#4c72b0",  # Standard academic blue
    edgecolor="black",
)

# 4. Labels and Title
plt.title(
    "Distribution of Observed Penetrance Across TTR Variants", fontsize=14, pad=15
)
plt.xlabel("Calculated Penetrance (%)", fontsize=12)
plt.ylabel("Number of Unique Variants", fontsize=12)

# Optional: Add a text annotation for the count of 100% penetrance variants
high_penetrance_count = data_to_plot[data_to_plot >= 99].count()
plt.text(
    95,
    5,
    f"High Penetrance\n(>99%): {high_penetrance_count}",
    ha="right",
    va="center",
    fontsize=10,
    color="darkred",
)

plt.tight_layout()
plt.show()
