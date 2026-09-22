# Lights Out — resolución mediante sistemas lineales sobre Z₂

Implementación en Python del punto 6 del proyecto de Álgebra Aplicada: una
función que recibe el estado inicial de un tablero de *Lights Out* de tamaño
n × n y devuelve el vector de pulsaciones que apaga todas las luces.

**Universidad Católica del Uruguay** — Facultad de Ingeniería y Tecnologías
Ingeniería en Informática — Álgebra Aplicada

Martina Caballero · Valentín Curbelo · Santiago Umpiérrez
Docentes: Catarina García y María José Gómez

---

## Qué hace

El juego se modela como el sistema lineal **A x = b** sobre el cuerpo Z₂, donde:

- **x** es el vector incógnita de tamaño n²: `x[i] = 1` significa que hay que
  presionar la luz *i*.
- **b** es el estado inicial del tablero leído por filas.
- **A** es la matriz de vecindades: `A[i][j] = 1` si presionar la luz *j*
  cambia el estado de la luz *i*.

El sistema se resuelve con una escalerización gaussiana adaptada a la
aritmética módulo 2, en la que la única transformación elemental empleada es
`Fi → Fi + Fj`, tal como pide la consigna.

## Requisitos

Solo Python 3.9 o superior. **No usa ninguna librería externa**: ni NumPy, ni
SymPy, ni nada que resuelva el sistema por nosotros. La escalerización está
implementada a mano.

## Uso desde la terminal

```bash
python lights_out.py
```

El programa pide el tamaño del tablero y luego cada fila. Si se presiona Enter
sin escribir nada, usa el tablero de ejemplo de la consigna.

Ejemplo de sesión completa:

```
==============================================================
Lights Out: resolución por sistemas lineales sobre Z_2
==============================================================
Tamaño del tablero n (Enter para usar el ejemplo 3x3): 3

Ingrese las 3 filas del tablero (1 = encendida, 0 = apagada).
Por ejemplo, una fila válida sería: 0 0 0

  Fila 1: 0 1 0
  Fila 2: 1 1 0
  Fila 3: 0 0 1

Tablero inicial
0 1 0
1 1 0
0 0 1

Solución (vector)
[0, 0, 0, 0, 1, 0, 0, 0, 1]

Luces que hay que presionar (1 = presionar)
0 0 0
0 1 0
0 0 1

Tablero tras aplicar esas pulsaciones
0 0 0
0 0 0
0 0 0

¿Quedaron todas las luces apagadas? sí
```

Las filas se pueden escribir con espacios (`0 1 0`) o sin ellos (`010`).

## Uso como módulo

La función que pide la consigna es `resolver_lights_out`:

```python
from lights_out import resolver_lights_out, verificar

tablero = [[0, 1, 0],
           [1, 1, 0],
           [0, 0, 1]]

x = resolver_lights_out(tablero)   # [0, 0, 0, 0, 1, 0, 0, 0, 1]
verificar(tablero, x)              # True
```

Devuelve `None` cuando el tablero no tiene solución. Por ejemplo, el tablero
de 5 × 5 que la consigna usa para ilustrar las reglas:

```python
tablero = [[0, 0, 1, 1, 1],
           [0, 1, 1, 1, 0],
           [1, 0, 0, 1, 1],
           [1, 1, 0, 0, 1],
           [0, 1, 1, 1, 0]]

resolver_lights_out(tablero)       # None: el sistema es inconsistente
```

## Funciones principales

| Función | Qué hace |
|---|---|
| `resolver_lights_out(tablero)` | **La función pedida.** Devuelve el vector solución o `None` |
| `construir_matriz_sistema(n)` | Construye la matriz A de vecindades, de n² × n² |
| `vector_terminos_independientes(t)` | Convierte el tablero en el vector b |
| `escalerizar(M, columnas)` | Escalerización gaussiana en Z₂ (solo `Fi → Fi + Fj`) |
| `resolver_sistema(A, b)` | Devuelve una solución particular y una base del espacio nulo |
| `todas_las_soluciones(tablero)` | Devuelve el conjunto solución completo |
| `presionar(tablero, fila, col)` | Simula una pulsación (para verificar) |
| `verificar(tablero, x)` | Comprueba que el vector realmente apaga el tablero |

Las funciones de simulación (`presionar`, `aplicar_solucion`, `verificar`) son
independientes del resolvedor: permiten comprobar los resultados jugando de
verdad, sin usar álgebra lineal.

## Cómo está organizado el archivo

`lights_out.py` está dividido en siete secciones comentadas:

1. Correspondencia entre las luces del tablero y las incógnitas del sistema
2. Construcción del sistema lineal A x = b
3. Escalerización gaussiana en Z₂
4. Resolución del sistema
5. Simulación del juego, para verificar las soluciones
6. Utilidades de presentación
7. Uso interactivo desde la terminal

## Nota sobre la resolubilidad

No todos los tableros tienen solución. Eso depende del rango de A, que a su vez
depende solo de n:

| n | rango A | dim Nul A | Soluciones por tablero | Tableros resolubles |
|---|---|---|---|---|
| 3 | 9 | 0 | 1 | todos |
| 4 | 12 | 4 | 16 | 1 de cada 16 |
| 5 | 23 | 2 | 4 | 1 de cada 4 |
| 6 | 36 | 0 | 1 | todos |

Cuando A es invertible, todo tablero se gana de una única manera. Cuando no lo
es, solo una fracción de los tableros tiene solución, y cada uno admite varias.
