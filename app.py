import math
import re

import numpy as np
import sympy as sp
from flask import Flask, render_template, request
from sympy.parsing.sympy_parser import (
    convert_xor,
    implicit_multiplication_application,
    parse_expr,
    standard_transformations,
)

app = Flask(__name__)

x = sp.Symbol("x", real=True)
TRANSFORMACIONES = standard_transformations + (
    implicit_multiplication_application,
    convert_xor,
)

# Lista blanca: lo único que se puede escribir dentro de una función
PERMITIDOS = {
    "x": x, "pi": sp.pi, "e": sp.E,
    "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
    "asin": sp.asin, "acos": sp.acos, "atan": sp.atan,
    "sinh": sp.sinh, "cosh": sp.cosh, "tanh": sp.tanh,
    "exp": sp.exp, "log": sp.log, "ln": sp.log,
    "sqrt": sp.sqrt, "abs": sp.Abs,
}
EJEMPLOS = ["x**2 - 4", "x^3 - 3x", "sin(x)", "exp(-x^2)", "1/x", "tan(x)"]
PUNTOS = 400


def interpretar(texto):
    """Convierte el texto en una expresión de sympy, validando antes."""
    if not texto:
        raise ValueError("Escribí una función, por ejemplo x**2 - 4.")
    if len(texto) > 100:
        raise ValueError("La función es demasiado larga.")
    if not re.fullmatch(r"[0-9A-Za-z+\-*/^().\s]+", texto):
        raise ValueError("Hay caracteres que no están permitidos.")
    for nombre in re.findall(r"[A-Za-z]\w*", texto):
        if nombre not in PERMITIDOS:
            raise ValueError(
                f"No conozco «{nombre}». Usá x y funciones como sin, cos, exp, log, sqrt."
            )
    try:
        expr = parse_expr(texto, local_dict=PERMITIDOS, transformations=TRANSFORMACIONES)
    except Exception:
        raise ValueError("No pude entender esa expresión. Revisá paréntesis y operadores.")
    if expr.free_symbols - {x}:
        raise ValueError("Solo se puede usar la variable x.")
    return expr


def crear_evaluador(expr):
    """Devuelve una función que evalúa expr en un número (o None si no se puede)."""
    f = sp.lambdify(x, expr, modules="mpmath")

    def evaluar(v):
        try:
            y = float(f(v))
            return y if math.isfinite(y) else None
        except Exception:
            return None

    return evaluar


def raices_numericas(evaluar, xs, ys):
    """Busca cambios de signo y los afina por bisección."""
    encontradas = []
    for i in range(len(xs) - 1):
        y0, y1 = ys[i], ys[i + 1]
        if y0 is None or y1 is None:
            continue
        if y0 == 0:
            encontradas.append(xs[i])
            continue
        if y0 * y1 < 0:
            a, b, fa = xs[i], xs[i + 1], y0
            for _ in range(40):
                m = (a + b) / 2
                fm = evaluar(m)
                if fm is None:
                    break
                if fa * fm <= 0:
                    b = m
                else:
                    a, fa = m, fm
            r = (a + b) / 2
            v = evaluar(r)
            if v is not None and abs(v) < 1e-6:  # descarta asíntotas
                encontradas.append(r)
    return encontradas[:30]


def calcular_raices(expr, evaluar, xs, ys):
    """Devuelve (textos en LaTeX, valores numéricos, si son exactas)."""
    try:
        sol = sp.solveset(sp.Eq(expr, 0), x, domain=sp.S.Reals)
    except Exception:
        sol = None

    if sol == sp.S.EmptySet:
        return [], [], True
    if sol == sp.S.Reals:
        return [r"\in \mathbb{R}"], [], True
    if isinstance(sol, sp.FiniteSet):
        raices = sorted(sol, key=lambda r: float(r))
        textos = []
        for r in raices:
            t = sp.latex(r)
            if not r.is_Rational:
                t += rf" \approx {float(r):.4f}"
            textos.append(t)
        return textos, [float(r) for r in raices], True

    # Si sympy no da una lista finita (ej. sin(x)), buscamos numéricamente
    aprox = raices_numericas(evaluar, xs, ys)
    return [rf"\approx {r:.4f}" for r in aprox], aprox, False


def serie(xs, ys):
    return [
        {"x": round(a, 4), "y": None if b is None else round(b, 6)}
        for a, b in zip(xs, ys)
    ]


def rango_y(*series):
    """Rango vertical razonable: ignora asíntotas pero respeta los picos reales."""
    vals = [v for s in series if s for v in s if v is not None]
    if not vals:
        return -1.0, 1.0
    lo, hi = np.percentile(vals, [2, 98])
    span = max(hi - lo, 1e-9)
    # Si el mínimo/máximo real está cerca del grueso de los datos, es un pico de verdad
    real_lo, real_hi = min(vals), max(vals)
    if real_lo >= lo - 0.5 * span:
        lo = real_lo
    if real_hi <= hi + 0.5 * span:
        hi = real_hi
    if hi - lo < 1e-9:
        lo, hi = lo - 1, hi + 1
    pad = (hi - lo) * 0.15
    return float(lo - pad), float(hi + pad)


def resolver(texto, xmin_txt, xmax_txt):
    try:
        xmin, xmax = float(xmin_txt), float(xmax_txt)
    except ValueError:
        raise ValueError("Los límites del eje x tienen que ser números.")
    if not xmin < xmax or max(abs(xmin), abs(xmax)) > 1000:
        raise ValueError("Poné límites válidos: mínimo menor al máximo, dentro de ±1000.")

    f = interpretar(texto)
    d = sp.diff(f, x)
    F = sp.integrate(f, x)
    integral_ok = not F.has(sp.Integral)

    xs = np.linspace(xmin, xmax, PUNTOS).tolist()
    ev_f = crear_evaluador(f)
    ev_d = crear_evaluador(d)
    ys_f = [ev_f(v) for v in xs]
    ys_d = [ev_d(v) for v in xs]
    if all(y is None for y in ys_f):
        raise ValueError("No pude graficar esa función en ese rango.")

    ys_i = None
    if integral_ok:
        ev_i = crear_evaluador(F)
        ys_i = [ev_i(v) for v in xs]

    textos, valores, exactas = calcular_raices(f, ev_f, xs, ys_f)
    en_rango = [r for r in valores if xmin <= r <= xmax]
    ymin, ymax = rango_y(ys_f, ys_d)

    return {
        "f_tex": sp.latex(f),
        "d_tex": sp.latex(d),
        "i_tex": sp.latex(F) + " + C" if integral_ok else None,
        "raices": textos,
        "raices_exactas": exactas,
        "grafico": {
            "f": serie(xs, ys_f),
            "d": serie(xs, ys_d),
            "i": serie(xs, ys_i) if ys_i else None,
            "raices": [{"x": round(r, 6), "y": 0} for r in en_rango],
            "ymin": ymin,
            "ymax": ymax,
        },
    }


@app.route("/")
def index():
    texto = request.args.get("f", "").strip()
    xmin = request.args.get("xmin", "-10")
    xmax = request.args.get("xmax", "10")
    contexto = {
        "texto": texto, "xmin": xmin, "xmax": xmax,
        "ejemplos": EJEMPLOS, "error": None, "resultado": None,
    }
    if texto:
        try:
            contexto["resultado"] = resolver(texto, xmin, xmax)
        except ValueError as e:
            contexto["error"] = str(e)
        except Exception:
            contexto["error"] = "Algo salió mal al calcular. Probá con otra función."
    return render_template("index.html", **contexto)


if __name__ == "__main__":
    app.run(debug=True)