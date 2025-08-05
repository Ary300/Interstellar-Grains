import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_results(results_file="results/science_campaign_v1.csv", output_dir="results/plots"):
    df = pd.read_csv(results_file)
    os.makedirs(output_dir, exist_ok=True)
    has_agg = "total_h2_formed_mean" in df.columns
    def col(base):
        return f"{base}_mean" if has_agg and f"{base}_mean" in df.columns else base
    def err(base):
        return f"{base}_ci95" if has_agg and f"{base}_ci95" in df.columns else None
    plt.figure(figsize=(10, 6))
    for density in sorted(df["h_gas_density_cm3"].unique()):
        subset = df[df["h_gas_density_cm3"] == density].sort_values("surface_temperature_k")
        y = subset[col("total_h2_formed")]
        yerr = subset[err("total_h2_formed")] if err("total_h2_formed") else None
        plt.errorbar(subset["surface_temperature_k"], y, yerr=yerr, marker="o", capsize=3, label=f"H Gas Density: {density}")
    plt.xlabel("Surface Temperature (K)")
    plt.ylabel("Total H2 Formed" + (" (mean ± 95% CI)" if has_agg else ""))
    plt.title("Total H2 Formed vs. Surface Temperature")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "h2_formed_vs_surface_temp.png"))
    plt.close()
    plt.figure(figsize=(10, 6))
    for temp in sorted(df["surface_temperature_k"].unique()):
        subset = df[df["surface_temperature_k"] == temp].sort_values("h_gas_density_cm3")
        y = subset[col("total_h2_formed")]
        yerr = subset[err("total_h2_formed")] if err("total_h2_formed") else None
        plt.errorbar(subset["h_gas_density_cm3"], y, yerr=yerr, marker="o", capsize=3, label=f"Surface Temp: {temp} K")
    plt.xlabel("H Gas Density (cm^-3)")
    plt.ylabel("Total H2 Formed" + (" (mean ± 95% CI)" if has_agg else ""))
    plt.title("Total H2 Formed vs. H Gas Density")
    plt.xscale("log")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "h2_formed_vs_h_gas_density.png"))
    plt.close()
    mech_cols = ["h2_formed_LH", "h2_formed_ER", "h2_formed_UV"]
    value_cols = [col(m) for m in mech_cols]
    plot_df = df.melt(id_vars=["surface_temperature_k", "h_gas_density_cm3"], value_vars=value_cols, var_name="Mechanism", value_name="Value")
    plot_df["Mechanism"] = plot_df["Mechanism"].str.replace("_mean$", "", regex=True)
    import seaborn as sns
    g = sns.catplot(data=plot_df, x="surface_temperature_k", y="Value", hue="Mechanism", col="h_gas_density_cm3", kind="bar", palette="viridis", errorbar=None)
    g.fig.suptitle("Contribution of H2 Formation Mechanisms", y=1.02)
    g.tight_layout()
    g.savefig(os.path.join(output_dir, "h2_formation_mechanisms_contribution.png"))
    plt.close()
    print(f"Plots saved to {output_dir}")

if __name__ == "__main__":
    plot_results()
