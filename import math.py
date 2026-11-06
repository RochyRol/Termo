import math

R = 8.314
N = 1.0
CV = 2.5 * R
CP = CV + R
GAMMA = CV / CP


def isobaric_work(P, V1, V2):
    return P * (V1 - V2)


def isothermal_work(n, T, V1, V2):
    return n * R * T * math.log10(V2 / V1)


def adiabatic_temperature(T1, V1, V2):
    return T1 * (V1 / V2) ** GAMMA


def isobaric_process(n, P, V1, V2):
    T1 = P * V1 / (n * R)
    T2 = P * V1 / (n * R)
    W = isobaric_work(P, V1, V2)
    dU = n * CV * (T2 - T1)
    Q = dU + W
    return Q, W, dU


def isochoric_process(n, T1, T2):
    W = 0.0
    dU = n * CV * (T2 - T1)
    Q = dU + W
    return Q, W, dU


def isothermal_process(n, T, V1, V2):
    W = isothermal_work(n, T, V1, V2)
    dU = 0.0
    Q = dU + W
    return Q, W, dU


def adiabatic_process(n, T1, V1, V2):
    T2 = adiabatic_temperature(T1, V1, V2)
    Q = 0.0
    dU = n * CV * (T2 - T1)
    W = Q - dU
    return Q, W, dU


def print_results(name, Q, W, dU):
    residual = Q - W - dU
    print(f"{name:<12} Q = {Q:9.2f} J   W = {W:9.2f} J   dU = {dU:9.2f} J   residual = {residual:.2e}")


def main():
    P = 101325.0
    V1 = 0.024
    V2 = 0.048
    T0 = 300.0

    print_results("Isobaric", *isobaric_process(N, P, V1, V2))
    print_results("Isochoric", *isochoric_process(N, T0, 2 * T0))
    print_results("Isothermal", *isothermal_process(N, T0, V1, V2))
    print_results("Adiabatic", *adiabatic_process(N, T0, V1, V2))


if __name__ == "__main__":
    main()