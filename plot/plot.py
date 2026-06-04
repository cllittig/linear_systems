"""
Gera gráficos de escalabilidade a partir do CSV do benchmark.

Uso simples (um benchmark):
    python plot/plot.py plot/data.csv plot/data.png

Uso comparativo (múltiplas variantes):
    make comparar
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

COLORS = {
    "lu":              "#e15759",
    "cholesky":        "#f28e2b",
    "cg":              "#4e79a7",
    "gauss_seidel":    "#76b7b2",
    "lapack_lu":       "#59a14f",
    "lapack_cholesky": "#499894",
}

LABELS = {
    "lu":              "LU (nossa impl.)",
    "cholesky":        "Cholesky (nossa impl.)",
    "cg":              "Gradiente Conjugado",
    "gauss_seidel":    "Gauss-Seidel",
    "lapack_lu":       "LAPACK dgesv",
    "lapack_cholesky": "LAPACK dpotrf/dpotrs",
}

_LINESTYLES = ["-", "--", ":", "-."]
_MARKERS    = ["o", "s", "^", "D"]


def _estilo_por_variante(variantes):
    estilos = {}
    for i, v in enumerate(sorted(variantes)):
        estilos[v] = {
            "linestyle": _LINESTYLES[i % len(_LINESTYLES)],
            "marker":    _MARKERS[i % len(_MARKERS)],
        }
    return estilos


def load(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    return df


def _grupos(df: pd.DataFrame):
    """Iterável de (label, cor, estilo, group) — agnóstico à presença de variante."""
    tem_variante = "variante" in df.columns and df["variante"].nunique() > 1
    if tem_variante:
        estilos = _estilo_por_variante(df["variante"].unique())
        for (method, variante), group in df.groupby(["metodo", "variante"]):
            label = f"{LABELS.get(method, method)} ({variante})"
            yield label, COLORS.get(method), estilos[variante], group.sort_values("n")
    else:
        for method, group in df.groupby("metodo"):
            estilo = {"linestyle": "-", "marker": "o"}
            yield LABELS.get(method, method), COLORS.get(method), estilo, group.sort_values("n")


def plot_time(df: pd.DataFrame, ax: plt.Axes) -> None:
    for label, cor, estilo, group in _grupos(df):
        ax.plot(group["n"], group["tempo_ms"],
                color=cor, label=label, linewidth=2, **estilo)
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
    for label, cor, estilo, group in _grupos(df):
        ax.plot(group["n"], group["residuo"],
                color=cor, label=label, linewidth=2, **estilo)
    ax.set_yscale("log")
    ax.set_xlabel("n")
    ax.set_ylabel("‖Ax − b‖")
    ax.set_title("Resíduo × n")
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.legend()
    ax.grid(True, which="both", linestyle="--", alpha=0.4)


def main() -> None:
    if len(sys.argv) < 3:
        print("Uso: python plot.py <data.csv> <saida.png> [--split]")
        sys.exit(1)

    df = load(sys.argv[1])
    out = sys.argv[2]
    split = "--split" in sys.argv

    tem_variante = "variante" in df.columns and df["variante"].nunique() > 1
    titulo = ("Comparação de variantes" if tem_variante
              else "Comparação de métodos para sistemas lineares densos")

    if split:
        base = out.rsplit(".", 1)[0]
        ext  = out.rsplit(".", 1)[-1]

        fig_t, ax_t = plt.subplots(figsize=(7, 5))
        fig_t.suptitle(titulo, fontsize=13)
        plot_time(df, ax_t)
        fig_t.tight_layout()
        out_t = f"{base}_tempo.{ext}"
        fig_t.savefig(out_t, dpi=150)
        print(f"Salvo em {out_t}")

        fig_r, ax_r = plt.subplots(figsize=(7, 5))
        fig_r.suptitle(titulo, fontsize=13)
        plot_residual(df, ax_r)
        fig_r.tight_layout()
        out_r = f"{base}_residuo.{ext}"
        fig_r.savefig(out_r, dpi=150)
        print(f"Salvo em {out_r}")
    else:
        fig, (ax_time, ax_res) = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(titulo, fontsize=13)
        plot_time(df, ax_time)
        plot_residual(df, ax_res)
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        print(f"Salvo em {out}")


if __name__ == "__main__":
    main()
