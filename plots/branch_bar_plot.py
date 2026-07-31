import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'

selected_models = [
    "gpt-5-mini",
    "gemini-3.1-pro-preview",
    "deepseek-reasoner",
    "meta-llama-3.1-70b",
]

display_names = {
    "gpt-5-mini":             "GPT-5 Mini",
    "gemini-3.1-pro-preview": "Gemini 3.1 Pro",
    "deepseek-reasoner":      "DeepSeek R1",
    "meta-llama-3.1-70b":    "Llama 3.1 70B",
}

data = {
    "gpt-5-mini":             {"Chemical": 60.00, "Electrical": 64.89, "Mechanical": 66.22},
    "gemini-3.1-pro-preview": {"Chemical": 54.22, "Electrical": 70.67, "Mechanical": 71.33},
    "deepseek-reasoner":      {"Chemical": 64.89, "Electrical": 57.33, "Mechanical": 60.89},
    "meta-llama-3.1-70b":    {"Chemical": 38.00, "Electrical": 48.89, "Mechanical": 35.56},
}

branch_order = ["Chemical", "Electrical", "Mechanical"]

color_palettes = {
    "gpt-5-mini":             ["#6ab0e8", "#9ccbf0", "#c8e4f8"],
    "gemini-3.1-pro-preview": ["#72c47a", "#a3d9a5", "#cceecf"],
    "deepseek-reasoner":      ["#b09cd6", "#ccc0e6", "#e3dcf3"],
    "meta-llama-3.1-70b":    ["#e08060", "#eda888", "#f5cdb8"],
}

hatches = {"Chemical": "", "Electrical": "///", "Mechanical": "..."}

fig, ax = plt.subplots(figsize=(14, 8))

x = np.arange(len(selected_models))
width = 0.22
offsets = {"Chemical": -width, "Electrical": 0, "Mechanical": width}

for branch in branch_order:
    bar_colors = [color_palettes[m][branch_order.index(branch)] for m in selected_models]
    bars = ax.bar(
        x + offsets[branch],
        [data[m][branch] for m in selected_models],
        width,
        color=bar_colors,
        hatch=hatches[branch],
        edgecolor="#888888",
        linewidth=0.7,
        alpha=1.0,
    )
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            h + 1.2,
            f"{h:.1f}",
            ha="center", va="bottom",
            fontsize=19, fontweight="bold",
            color="black", fontname="Times New Roman",
        )

ax.set_xticks(x)
ax.set_xticklabels(
    [display_names[m] for m in selected_models],
    fontsize=27, fontweight="bold", fontname="Times New Roman",
)
ax.set_ylabel("Final Answer Accuracy (%)", fontsize=20, fontweight="bold", fontname="Times New Roman")
ax.set_ylim(0, 100)
ax.tick_params(axis="y", labelsize=16, colors="black")
ax.tick_params(axis="x", colors="black", pad=12)

ax.grid(axis="y", linestyle="--", alpha=0.4, linewidth=0.8, color="gray", zorder=0)
ax.grid(axis="x", alpha=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("black");   ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_color("black"); ax.spines["bottom"].set_linewidth(1.2)

legend_elements = [
    mpatches.Patch(facecolor="#cccccc", edgecolor="black", linewidth=1.2,
                   hatch=hatches["Chemical"],   label="Chemical Eng."),
    mpatches.Patch(facecolor="#cccccc", edgecolor="black", linewidth=1.2,
                   hatch=hatches["Electrical"], label="Electrical Eng."),
    mpatches.Patch(facecolor="#cccccc", edgecolor="black", linewidth=1.2,
                   hatch=hatches["Mechanical"], label="Mechanical Eng."),
]
legend = ax.legend(
    handles=legend_elements,
    title="Engineering Branch",
    loc="upper right",
    frameon=True, fancybox=False, shadow=False,
    borderpad=1.0, edgecolor="black", framealpha=1, facecolor="white",
    fontsize=15, title_fontsize=17,
    prop={"family": "Times New Roman", "size": 15},
)
legend.get_title().set_fontweight("bold")
legend.get_title().set_fontname("Times New Roman")
legend.get_frame().set_linewidth(1.5)

plt.tight_layout()
plt.savefig("Branch_Performance_Bar.pdf", dpi=300, bbox_inches="tight")
plt.savefig("Branch_Performance_Bar.png", dpi=300, bbox_inches="tight")
plt.show()
