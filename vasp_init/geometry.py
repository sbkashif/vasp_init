import math
import re
from typing import List


def _normalize_element(sym: str) -> str:
    s = re.sub(r"[^A-Za-z]", "", sym.strip())
    if not s:
        return "X"
    if len(s) == 1:
        return s.upper()
    return s[0].upper() + s[1:].lower()


def det3(m: List[List[float]]) -> float:
    (a,b,c),(d,e,f),(g,h,i) = m
    return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)


def mat_inv3(m: List[List[float]]) -> List[List[float]]:
    (a,b,c),(d,e,f),(g,h,i) = m
    A = e*i - f*h
    B = -(d*i - f*g)
    C = d*h - e*g
    D = -(b*i - c*h)
    E = a*i - c*g
    F = -(a*h - b*g)
    G = b*f - c*e
    H = -(a*f - c*d)
    I = a*e - b*d
    det = a*A + b*B + c*C
    if abs(det) < 1e-14:
        raise ValueError("Singular lattice matrix (det ~ 0)")
    inv = [
        [A/det, D/det, G/det],
        [B/det, E/det, H/det],
        [C/det, F/det, I/det],
    ]
    return inv


def mat_vec(m: List[List[float]], v: List[float]) -> List[float]:
    return [m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
            m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
            m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]]


def vec_mod1(frac: List[float]) -> List[float]:
    return [x - math.floor(x) for x in frac]
