"""
Gera gráficos de escalabilidade a partir do CSV do benchmark.

Uso:
    make benchmark > plot/data.csv
    cd plot && .venv/bin/python plot.py data.csv
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

COLORS = {
    "lu":           "#e15759",
    "cholesky":     "#f28e2b",
    "cg":           "#4e79a7",
    "gauss_seidel": "#76b7b2",
}

LABELS = {
    "lu":           "LU",
    "cholesky":     "Cholesky",
    "cg":           "Gradiente Conjugado",
    "gauss_seidel": "Gauss-Seidel",
}


def load(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    return df


def plot_time(df: pd.DataFrame, ax: plt.Axes) -> None:
    for method, group in df.groupby("metodo"):
        group = group.sort_values("n")
        ax.plot(
            group["n"], group["tempo_ms"],
            marker="o", label=LABELS.get(method, method),
            color=COLORS.get(method), linewidth=2,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("n")
    ax.set_ylabel("Tempo (ms)")
    ax.set_title("Tempo × n  (escala log-log)")
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.4)


def plot_residual(df: pd.DataFrame, ax: plt.Axes) -> None:
    for method, group in df.groupby("metodo"):
        group = group.sort_values("n")
        ax.plot(
            group["n"], group["residuo"],
            marker="s", label=LABELS.get(method, method),
            color=COLORS.get(method), linewidth=2,
        )
    ax.set_yscale("log")
    ax.set_xlabel("n")
    ax.set_ylabel("‖Ax − b‖")
    ax.set_title("Resíduo × n")
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.4)


def main() -> None:
    if len(sys.argv) < 2:
        print("Uso: python plot.py <data.csv>")
        sys.exit(1)

    df = load(sys.argv[1])

    fig, (ax_time, ax_res) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Comparação de métodos para sistemas lineares densos", fontsize=13)

    plot_time(df, ax_time)
    plot_residual(df, ax_res)

    fig.tight_layout()
    out = sys.argv[1].replace(".csv", ".png")
    fig.savefig(out, dpi=150)
    print(f"Salvo em {out}")


if __name__ == "__main__":
    main()
