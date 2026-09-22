"""
lights_out.py
-------------
Modelo algebraico y resolución del juego Lights Out sobre el cuerpo Z_2.

El tablero de n x n luces se modela como un sistema lineal A x = b con
n^2 ecuaciones y n^2 incógnitas, donde toda la aritmética se realiza
módulo 2 (1 + 1 = 0).

La resolución se hace con un algoritmo de escalerización gaussiana
adaptado a Z_2: la única transformación elemental de filas utilizada es
F_i -> F_i + F_j, ya que en Z_2 el escalamiento de filas es trivial
(el único escalar no nulo es 1) y el intercambio de filas puede
sustituirse por una suma.

Proyecto de Álgebra Aplicada - Universidad Católica del Uruguay.
"""

from __future__ import annotations

Tablero = list[list[int]]
Vector = list[int]
Matriz = list[list[int]]


# ---------------------------------------------------------------------------
# 1. Correspondencia entre las luces del tablero y las incógnitas del sistema
# ---------------------------------------------------------------------------

def indice_de_luz(fila: int, columna: int, n: int) -> int:
    """Devuelve el índice i (base 1) de la incógnita asociada a la luz a_{fila,columna}.

    Se usa la regla del enunciado: si i - 1 = n * q + r con 0 <= r < n,
    entonces la incógnita x_i corresponde a la luz a_{(q+1)(r+1)}.
    Los argumentos `fila` y `columna` se indican en base 1.
    """
    return n * (fila - 1) + (columna - 1) + 1


def luz_de_indice(i: int, n: int) -> tuple[int, int]:
    """Operación inversa de `indice_de_luz`: devuelve (fila, columna) en base 1."""
    q, r = divmod(i - 1, n)
    return q + 1, r + 1


def vecindad(fila: int, columna: int, n: int) -> list[tuple[int, int]]:
    """Luces cuyo estado cambia al presionar a_{fila,columna}.

    Son la propia luz presionada y sus vecinas ortogonales que caen
    dentro del tablero (arriba, abajo, izquierda y derecha).
    """
    candidatas = [
        (fila, columna),
        (fila - 1, columna),
        (fila + 1, columna),
        (fila, columna - 1),
        (fila, columna + 1),
    ]
    return [(f, c) for f, c in candidatas if 1 <= f <= n and 1 <= c <= n]


# ---------------------------------------------------------------------------
# 2. Construcción del sistema lineal A x = b
# ---------------------------------------------------------------------------

def construir_matriz_sistema(n: int) -> Matriz:
    """Construye la matriz de coeficientes A del sistema, de tamaño n^2 x n^2.

    La entrada A[i][j] vale 1 si al presionar la luz asociada a la incógnita
    x_{j+1} cambia el estado de la luz asociada a la ecuación i+1, y vale 0
    en caso contrario. La matriz resulta simétrica, porque la relación
    "ser vecinas" es simétrica.
    """
    N = n * n
    A = [[0] * N for _ in range(N)]
    for fila in range(1, n + 1):
        for columna in range(1, n + 1):
            i = indice_de_luz(fila, columna, n)
            for f_vecina, c_vecina in vecindad(fila, columna, n):
                j = indice_de_luz(f_vecina, c_vecina, n)
                A[i - 1][j - 1] = 1
    return A


def vector_terminos_independientes(tablero: Tablero) -> Vector:
    """Devuelve el vector b del sistema: el estado inicial del tablero leído por filas.

    La componente b_i vale 1 si la luz asociada a la ecuación i está
    encendida, es decir, si se le debe provocar un cambio de estado.
    """
    return [valor for fila in tablero for valor in fila]


def matriz_aumentada(A: Matriz, b: Vector) -> Matriz:
    """Devuelve la matriz aumentada [A | b]."""
    return [fila[:] + [termino] for fila, termino in zip(A, b)]


# ---------------------------------------------------------------------------
# 3. Escalerización gaussiana en Z_2
# ---------------------------------------------------------------------------

def sumar_filas(M: Matriz, destino: int, origen: int) -> None:
    """Aplica la transformación elemental F_destino -> F_destino + F_origen en Z_2.

    La suma en Z_2 coincide con el "o exclusivo" bit a bit, de modo que
    sumar dos veces la misma fila deja la matriz sin cambios.
    """
    fila_origen = M[origen]
    fila_destino = M[destino]
    for k in range(len(fila_destino)):
        fila_destino[k] = (fila_destino[k] + fila_origen[k]) % 2


def escalerizar(M: Matriz, columnas_incognitas: int) -> list[int]:
    """Lleva M a la forma escalonada reducida por filas trabajando en Z_2.

    `columnas_incognitas` indica cuántas columnas iniciales corresponden a
    incógnitas; las restantes se tratan como términos independientes.
    Devuelve la lista de columnas pivote encontradas (índices base 0).

    La matriz M se modifica en el lugar. La única transformación elemental
    empleada es F_i -> F_i + F_j: no hace falta escalar filas, porque el
    único coeficiente no nulo de Z_2 es 1, ni intercambiarlas, porque para
    instalar un pivote alcanza con sumarle a la fila actual alguna fila
    inferior que tenga un 1 en esa columna.
    """
    filas = len(M)
    columnas_pivote: list[int] = []
    fila_pivote = 0

    for columna in range(columnas_incognitas):
        if fila_pivote >= filas:
            break

        # Si no hay un 1 en la posición pivote, se lo instala sumando una
        # fila inferior que sí lo tenga. Si ninguna lo tiene, la columna no
        # es columna pivote y su incógnita será una variable libre.
        if M[fila_pivote][columna] == 0:
            candidata = next(
                (f for f in range(fila_pivote + 1, filas) if M[f][columna] == 1),
                None,
            )
            if candidata is None:
                continue
            sumar_filas(M, fila_pivote, candidata)

        # Se anulan el resto de las entradas de la columna, por debajo y por
        # encima del pivote, para obtener directamente la forma reducida.
        for fila in range(filas):
            if fila != fila_pivote and M[fila][columna] == 1:
                sumar_filas(M, fila, fila_pivote)

        columnas_pivote.append(columna)
        fila_pivote += 1

    return columnas_pivote


# ---------------------------------------------------------------------------
# 4. Resolución del sistema
# ---------------------------------------------------------------------------

def resolver_sistema(A: Matriz, b: Vector) -> tuple[Vector | None, list[Vector]]:
    """Resuelve A x = b en Z_2 por escalerización de la matriz aumentada.

    Devuelve el par (solución particular, base del espacio nulo de A):

    * si el sistema es inconsistente, la solución particular es None;
    * si es consistente, la solución particular se obtiene asignando el
      valor 0 a todas las variables libres, y la base del espacio nulo
      describe el conjunto solución completo como x = x_p + Nul A.
    """
    N = len(A)
    M = matriz_aumentada(A, b)
    columnas_pivote = escalerizar(M, N)

    # El sistema es inconsistente si aparece una fila de la forma
    # [0 0 ... 0 | 1], que corresponde a la ecuación 0 = 1.
    for fila in M:
        if fila[N] == 1 and all(valor == 0 for valor in fila[:N]):
            return None, []

    pivote_de_columna = {c: f for f, c in enumerate(columnas_pivote)}
    columnas_libres = [c for c in range(N) if c not in pivote_de_columna]

    # Solución particular: todas las variables libres valen 0 y cada
    # variable básica toma el término independiente de su fila.
    solucion = [0] * N
    for columna, fila in pivote_de_columna.items():
        solucion[columna] = M[fila][N]

    # Base del espacio nulo: se asigna el valor 1 a una variable libre por
    # vez y 0 a las demás, despejando las variables básicas.
    base_nulo: list[Vector] = []
    for libre in columnas_libres:
        generador = [0] * N
        generador[libre] = 1
        for columna, fila in pivote_de_columna.items():
            generador[columna] = M[fila][libre]
        base_nulo.append(generador)

    return solucion, base_nulo


def resolver_lights_out(tablero: Tablero) -> Vector | None:
    """Resuelve una partida de Lights Out.

    Recibe una matriz n x n de unos y ceros con el estado inicial del juego
    (1 = luz encendida, 0 = luz apagada) y devuelve un vector de n^2 unos y
    ceros, donde la componente x_i vale 1 si la luz correspondiente debe
    presionarse. Si el tablero no admite solución, devuelve None.
    """
    n = len(tablero)
    if any(len(fila) != n for fila in tablero):
        raise ValueError("El tablero debe ser una matriz cuadrada.")
    if any(valor not in (0, 1) for fila in tablero for valor in fila):
        raise ValueError("El tablero solo puede contener los valores 0 y 1.")

    A = construir_matriz_sistema(n)
    b = vector_terminos_independientes(tablero)
    solucion, _ = resolver_sistema(A, b)
    return solucion


def todas_las_soluciones(tablero: Tablero) -> list[Vector]:
    """Devuelve el conjunto solución completo del tablero dado.

    El conjunto solución es x_p + Nul A, de modo que si el espacio nulo
    tiene dimensión d el tablero admite exactamente 2^d soluciones.
    """
    n = len(tablero)
    A = construir_matriz_sistema(n)
    b = vector_terminos_independientes(tablero)
    particular, base_nulo = resolver_sistema(A, b)
    if particular is None:
        return []

    soluciones = [particular]
    for generador in base_nulo:
        soluciones += [
            [(u + v) % 2 for u, v in zip(sol, generador)] for sol in soluciones
        ]
    return soluciones


# ---------------------------------------------------------------------------
# 5. Simulación del juego (para verificar las soluciones encontradas)
# ---------------------------------------------------------------------------

def presionar(tablero: Tablero, fila: int, columna: int) -> Tablero:
    """Devuelve el tablero que resulta de presionar la luz a_{fila,columna}."""
    n = len(tablero)
    nuevo = [f[:] for f in tablero]
    for f, c in vecindad(fila, columna, n):
        nuevo[f - 1][c - 1] = (nuevo[f - 1][c - 1] + 1) % 2
    return nuevo


def aplicar_solucion(tablero: Tablero, x: Vector) -> Tablero:
    """Aplica al tablero todas las pulsaciones indicadas por el vector x."""
    n = len(tablero)
    resultado = [f[:] for f in tablero]
    for i in range(1, n * n + 1):
        if x[i - 1] == 1:
            fila, columna = luz_de_indice(i, n)
            resultado = presionar(resultado, fila, columna)
    return resultado


def esta_apagado(tablero: Tablero) -> bool:
    """Indica si todas las luces del tablero estan apagadas."""
    return all(valor == 0 for fila in tablero for valor in fila)


def verificar(tablero: Tablero, x: Vector) -> bool:
    """Comprueba que el vector x efectivamente gana la partida."""
    return esta_apagado(aplicar_solucion(tablero, x))


# ---------------------------------------------------------------------------
# 6. Utilidades de presentacion
# ---------------------------------------------------------------------------

def vector_a_tablero(x: Vector, n: int) -> Tablero:
    """Reacomoda un vector de tamaño n^2 como una matriz n x n leída por filas."""
    return [x[fila * n:(fila + 1) * n] for fila in range(n)]


def formatear_tablero(tablero: Tablero) -> str:
    """Devuelve una representación textual del tablero, una fila por línea."""
    return "\n".join(" ".join(str(valor) for valor in fila) for fila in tablero)


# ---------------------------------------------------------------------------
# 7. Uso interactivo desde la terminal
# ---------------------------------------------------------------------------

EJEMPLO_CONSIGNA: Tablero = [
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
]


def leer_tablero() -> Tablero:
    """Pide el tablero por teclado, fila por fila, y lo devuelve como matriz.

    Se ingresa primero el tamaño n y luego las n filas, cada una como una
    secuencia de n ceros y unos separados o no por espacios. Si el tamaño se
    deja en blanco, se usa el tablero de ejemplo de la consigna.
    """
    texto = input("Tamaño del tablero n (Enter para usar el ejemplo 3x3): ")
    if not texto.strip():
        return EJEMPLO_CONSIGNA

    n = int(texto)
    if n < 1:
        raise ValueError("El tamaño del tablero debe ser un entero positivo.")

    print(f"\nIngrese las {n} filas del tablero (1 = encendida, 0 = apagada).")
    print(f"Por ejemplo, una fila válida sería: {' '.join('0' * n)}\n")

    tablero: Tablero = []
    for fila in range(1, n + 1):
        while True:
            valores = [c for c in input(f"  Fila {fila}: ") if c in "01"]
            if len(valores) == n:
                tablero.append([int(c) for c in valores])
                break
            print(f"  Se esperaban {n} valores de 0 o 1. Intente de nuevo.")
    return tablero


def informar(tablero: Tablero) -> None:
    """Resuelve el tablero e imprime el resultado por pantalla."""
    n = len(tablero)

    print("\nTablero inicial")
    print(formatear_tablero(tablero))

    x = resolver_lights_out(tablero)
    if x is None:
        print("\nEste tablero NO tiene solución: el sistema es inconsistente,")
        print("de modo que ninguna combinación de pulsaciones lo apaga.")
        return

    print("\nSolución (vector)")
    print(x)

    print("\nLuces que hay que presionar (1 = presionar)")
    print(formatear_tablero(vector_a_tablero(x, n)))

    print("\nTablero tras aplicar esas pulsaciones")
    print(formatear_tablero(aplicar_solucion(tablero, x)))

    print("\n¿Quedaron todas las luces apagadas?",
          "sí" if verificar(tablero, x) else "no")

    cantidad = len(todas_las_soluciones(tablero))
    if cantidad > 1:
        print(f"\nEste tablero admite {cantidad} soluciones distintas; arriba se")
        print("muestra una de ellas.")


if __name__ == "__main__":
    print("=" * 62)
    print("Lights Out: resolución por sistemas lineales sobre Z_2")
    print("=" * 62)
    informar(leer_tablero())
