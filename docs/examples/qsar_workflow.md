# Example: QSAR Workflow

This example outlines a basic Quantitative Structure-Activity Relationship (QSAR) workflow using the results from `pyrene-dimer-analyzer`.

**Goal**: To build a model that predicts an experimental property (e.g., excimer-to-monomer emission ratio) based on the calculated geometric descriptors.

## 1. Generate Geometric Descriptors

First, run the analysis on all your conformer files to generate the geometric data.

```bash
pyrene-analyze analyze *.sdf -o geometric_descriptors.csv
```

## 2. Prepare Experimental Data

Next, you need a file with your experimental data. For this example, let's assume you have a file named `experimental_data.csv` with the following format:

```csv
molecule,excimer_ratio
Et,0.1
iPr,0.3
cHex,0.5
tBu,0.8
```

## 3. Build the QSAR Model

This Python script will load both datasets, merge them, and build a simple regression model.

```python
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

# Load the geometric descriptors
geom_df = pd.read_csv("geometric_descriptors.csv")

# Load the experimental data
exp_df = pd.read_csv("experimental_data.csv")

# Calculate the average geometric properties for each molecule
avg_geom_df = geom_df.groupby("molecule").mean().reset_index()

# Merge the geometric and experimental data
merged_df = pd.merge(avg_geom_df, exp_df, on="molecule")

# Define features (X) and target (y)
features = ["plane_angle_deg", "interplane_distance_A", "pi_overlap_pct"]
X = merged_df[features]
y = merged_df["excimer_ratio"]

# Scale the features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Build a linear regression model
model = LinearRegression()
model.fit(X_scaled, y)

# Print the model coefficients to see feature importance
print("QSAR Model Coefficients:")
for feature, coef in zip(features, model.coef_):
    print(f"  - {feature}: {coef:.3f}")

# You can now use this model to predict the excimer ratio for new molecules
# new_molecule_features = ...
# new_X_scaled = scaler.transform(new_molecule_features)
# predicted_ratio = model.predict(new_X_scaled)
```

This example provides a starting point for a more complex QSAR analysis. You can expand on this by using more advanced models (e.g., `RandomForestRegressor`), performing cross-validation, and engineering more complex features from the geometric data.
