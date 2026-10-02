# ∑ Graficador de funciones

App web en Flask que grafica una función y calcula su derivada, su integral y sus raíces con sympy.

## Qué hace

- Escribís una función (por ejemplo `x**2 - 4` o `sin(x)`) y elegís el rango del eje x.
- Muestra la derivada, la integral indefinida y las raíces con notación matemática (KaTeX).
- Grafica la función, su derivada y su integral con Chart.js.
- Valida el texto con una lista blanca de funciones permitidas, para no ejecutar código arbitrario.

## Tecnologías

Python · Flask · sympy · numpy · KaTeX · Chart.js

## Cómo correrlo

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python app.py
```

Después abrí http://127.0.0.1:5000 en el navegador.

## Limitaciones conocidas

- La escala vertical es automática: en funciones con asíntotas (como tan(x)) los picos de la derivada pueden aplastar el gráfico.
- Las ramas de funciones como 1/x se unen con una línea vertical a través de la asíntota.
- Sympy devuelve log(x) como integral de 1/x (lo correcto es ln|x|, válido para x ≠ 0).

## Capturas

![Resultados](capturas/resultados.png)
![Gráfico](capturas/grafico.png)