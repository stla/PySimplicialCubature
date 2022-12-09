from math import factorial, cos, acos, sqrt, pi
import numpy as np
from sys import exit
from sympy import Poly


def _SMPCHC(ND, NF, MXFS, EA, ER, SBS, KEY):
    FL = 0
    if (KEY < 0) | (KEY > 4):
        FL = 2
    if ND < 2:
        FL = 3
    if NF < 1:
        FL = 4
    if (EA < 0) & (ER < 0):
        FL = 5
    RULCLS = 0
    if FL == 0:
        if KEY == 0:
            RULCLS = (ND + 4) * (ND + 3) * (ND + 2) / 6 + (ND + 2) * (ND + 1)
        elif KEY == 1:
            RULCLS = 2 * ND + 3
        elif KEY == 2:
            RULCLS = (ND + 3) * (ND + 2) / 2 + 2 * (ND + 1)
        elif KEY == 3:
            RULCLS = (ND + 4) * (ND + 3) * (ND + 2) / 6 + (ND + 2) * (ND + 1)
        elif KEY == 4:
            RULCLS = (ND + 5) * (ND + 4) * (ND + 3) * (ND + 2) / 24 + 5 * (ND + 2) * (
                ND + 1
            ) / 2
        if MXFS < SBS * RULCLS:
            FL = 1
    return FL, RULCLS


def _SMPDFS(ND, NF, F, TOP, SBS, VRTS):
    SBS = SBS-1
    CUTTF = 2
    CUTTB = 8
    IS = 0
    JS = 1
    DFMX = 0
    EMX = 0
    V = VRTS[TOP]
    _, n = V.shape
    CN = np.apply_along_axis(np.sum, 1, V) / n
    FC = F(CN)
    DFMD = np.sum(np.abs(FC))
    FRTHDF = np.zeros((ND, ND + 1))
    for I in range(ND):
        for J in range(I + 1, ND + 1):
            H = 2 * (V[:, I] - V[:, J]) / (5 * (ND + 1))
            EWD = np.sum(np.abs(H))
            if EWD >= EMX:
                IE = I
                JE = J
                EMX = EWD
            DFR = np.sum(
                np.abs(
                    F(CN - 2 * H) + F(CN + 2 * H) + 6 * FC - 4 * (F(CN - H) + F(CN + H))
                )
            )
            if (DFMD + DFR / 8) == DFMD:
                DFR = 0
            DFR = DFR * EWD
            if DFR >= DFMX:
                IT = IS
                JT = JS
                DFNX = DFMX
                IS = I
                JS = J
                DFMX = DFR
            else:
                if DFR >= DFNX:  # ??? DFNX not defined !
                    IT = I
                    JT = J
                    DFNX = DFR
            FRTHDF[I, J] = DFR
    if DFNX > DFMX / CUTTF:
        NEW = 4
    else:
        NEW = 3
        if DFMX == 0:
            IS = IE
            JS = JE
        else:
            DFSMX = 0
            for L in range(ND + 1):
                if (L != IS) & (L != JS):
                    IT = np.min([L, IS, JS])
                    JT = np.max([L, IS, JS])
                    LT = IS + JS + L - IT - JT
                    DFR = FRTHDF[IT, LT] + FRTHDF[LT, JT]
                    if DFR >= DFSMX:
                        DFSMX = DFR
                        LS = L
            DIFIL = FRTHDF[min(IS, LS), max(IS, LS)]
            DIFLJ = FRTHDF[min(JS, LS), max(JS, LS)]
            DFNX = DIFIL + DIFLJ - min(DIFIL, DIFLJ)
            if (DFMX / CUTTB < DFNX) & (DIFIL > DIFLJ):
                IT = IS
                IS = JS
                JS = IT
    VV = np.repeat([V], NEW - 1, axis=0)
    VRTS = np.concatenate((VRTS, VV))
    VTI = V[:, IS]
    VTJ = V[:, JS]
    if NEW == 4:
        VRTS[TOP][:, JS] = (VTI + VTJ) / 2
        VRTS[SBS + 1][:, IS] = VTI
        VRTS[SBS + 1][:, JS] = (VTI + VTJ) / 2
        VRTS[SBS + 2][:, IS] = (VTI + VTJ) / 2
        VRTS[SBS + 2][:, JS] = VTJ
        VRTS[SBS + 3][:, IS] = (VTI + VTJ) / 2
        VRTS[SBS + 3][:, JS] = VTJ
        VTI = np.copy(VRTS[TOP][:, IT])
        VTJ = np.copy(VRTS[TOP][:, JT])
        VRTS[TOP][:, JT] = (VTI + VTJ) / 2
        VRTS[SBS + 1][:, IT] = (VTI + VTJ) / 2
        VRTS[SBS + 1][:, JT] = VTJ
        VTI = np.copy(VRTS[SBS + 2][:, IT])
        VTJ = np.copy(VRTS[SBS + 2][:, JT])
        VRTS[SBS + 2][:, JT] = (VTI + VTJ) / 2
        VRTS[SBS + 3][:, IT] = (VTI + VTJ) / 2
        VRTS[SBS + 3][:, JT] = VTJ
    else:
        VRTS[TOP][:, JS] = (2 * VTI + VTJ) / 3
        VRTS[SBS + 1][:, IS] = (2 * VTI + VTJ) / 3
        if DFMX / CUTTF < DFNX:
            VRTS[SBS + 1][:, JS] = VTJ
            VRTS[SBS + 2][:, IS] = (2 * VTI + VTJ) / 3
            VRTS[SBS + 2][:, JS] = VTJ
            VTJ = VRTS[SBS + 1][:, JS]
            VTL = VRTS[SBS + 1][:, LS]
            VRTS[SBS + 1][:, LS] = (VTJ + VTL) / 2
            VRTS[SBS + 2][:, JS] = (VTJ + VTL) / 2
            VRTS[SBS + 2][:, LS] = VTL
        else:
            VRTS[SBS + 1][:, JS] = (VTI + 2 * VTJ) / 3
            VRTS[SBS + 2][:, IS] = (VTI + 2 * VTJ) / 3
            VRTS[SBS + 2][:, JS] = VTJ
    return VRTS, NEW


def _SMPSMS(N, VERTEX, NF, F, G):
    SYMSMS = np.zeros(NF)
    G = np.flip(np.sort(G))
    pr = True
    while pr:
        SYMSMS = SYMSMS + F(VERTEX @ G)
        pr = False
        for I in range(1, N + 1):
            GI = G[I]
            if G[I - 1] > GI:
                IX = I - 1
                d, _ = divmod(IX + 2, 2)
                for L in range(d):  # !!!!!!!!!!
                    GL = G[L]
                    if GL <= GI:
                        IX = IX - 1
                    G[L] = G[I - L - 1]
                    G[I - L - 1] = GL
                    if G[L] > GI:
                        LX = L
                if G[IX] <= GI:
                    IX = LX
                G[I] = G[IX]
                G[IX] = GI
                pr = True
                break
    return np.reshape(SYMSMS, (1, 1))


def _SMPRMS(N, KEY):
    if KEY == 1:
        RLS = 3
        GMS = 2
        WTS = 3
    elif KEY == 2:
        RLS = 5
        GMS = 4
        WTS = 6
    elif KEY == 3:
        RLS = 7
        GMS = 7
        WTS = 11
    elif KEY == 4:
        RLS = 7
        GMS = 12
        WTS = 21
        if N == 2:
            GMS = 11
            WTS = 20
    W = np.zeros((WTS, RLS))
    PTS = np.zeros(WTS)
    G = np.zeros((N + 1, WTS))
    NP = N + 1
    N2 = NP * (N + 2)
    N4 = N2 * (N + 3) * (N + 4)
    N6 = N4 * (N + 5) * (N + 6)
    N8 = N6 * (N + 7) * (N + 8)
    G[:, 0] = 1 / NP
    PTS[0] = 1
    R1 = (N + 4 - sqrt(15)) / (N * N + 8 * N + 1)
    S1 = 1 - N * R1
    L1 = S1 - R1
    G[0, GMS] = S1
    G[1:NP, GMS] = R1
    PTS[GMS] = NP
    IW = RLS
    if KEY < 4:
        W[0, IW - 1] = 1
        IW = IW - 1
        W[GMS, IW - 1] = 1 / NP
        IW = IW - 1
    G[0, 1] = 3 / (N + 3)
    G[1:NP, 1] = 1 / (N + 3)
    PTS[1] = NP
    W[1, IW - 1] = (N + 3) ** 3 / (4 * N2 * (N + 3))
    if KEY > 1:
        IW = IW - 1
        if N == 2:
            L2 = 0.62054648267200632589046034361711
            L1 = -sqrt(1 / 2 - L2 ** 2)
            R1 = (1 - L1) / 3
            S1 = 1 - 2 * R1
            G[0, GMS] = S1
            G[1:NP, GMS] = R1
            PTS[GMS] = 3
            W[GMS, IW - 1] = 1 / 6
            R2 = (1 - L2) / 3
            S2 = 1 - 2 * R2
            G[0, GMS + 1] = S2
            G[1:NP, GMS + 1] = R2
            PTS[GMS + 1] = 3
            W[GMS + 1, IW - 1] = 1 / 6
        else:
            R2 = (N + 4 + sqrt(15)) / (N * N + 8 * N + 1)
            S2 = 1 - N * R2
            L2 = S2 - R2
            G[0, GMS + 1] = S2
            G[1:NP, GMS + 1] = R2
            PTS[GMS + 1] = NP
            W[GMS + 1, IW - 1] = (2 / (N + 3) - L1) / (N2 * (L2 - L1) * L2 ** 2)
            W[GMS, IW - 1] = (2 / (N + 3) - L2) / (N2 * (L1 - L2) * L1 ** 2)
        IW = IW - 1
        G[0, 2] = 5 / (N + 5)
        G[1:NP, 2] = 1 / (N + 5)
        PTS[2] = NP
        G[0, 3] = 3 / (N + 5)
        G[1, 3] = 3 / (N + 5)
        G[2:NP, 3] = 1 / (N + 5)
        PTS[3] = NP * N / 2
        W[1, IW - 1] = -((N + 3) ** 5) / (16 * N4)
        W[2:4, IW - 1] = (N + 5) ** 5 / (16 * N4 * (N + 5))
    if KEY > 2:
        IW = IW - 1
        U1 = (N + 7 + 2 * sqrt(15)) / (N * N + 14 * N - 11)
        V1 = (1 - (N - 1) * U1) / 2
        D1 = V1 - U1
        G[0, GMS + 2] = V1
        G[1, GMS + 2] = V1
        G[2:NP, GMS + 2] = U1
        PTS[GMS + 2] = ((N + 1) * N) / 2
        U2 = (N + 7 - 2 * sqrt(15)) / (N * N + 14 * N - 11)
        V2 = (1 - (N - 1) * U2) / 2
        D2 = V2 - U2
        G[0, GMS + 3] = V2
        G[1, GMS + 3] = V2
        G[2:NP, GMS + 3] = U2
        PTS[GMS + 3] = ((N + 1) * N) / 2
        if N == 2:
            W[GMS + 2, IW - 1] = (155 - sqrt(15)) / 1200
            W[GMS + 3, IW - 1] = (155 + sqrt(15)) / 1200
            W[0, IW - 1] = 1 - 3 * (W[GMS + 2, IW - 1] + W[GMS + 3, IW - 1])
        elif N == 3:
            W[GMS, IW - 1] = (2665 + 14 * sqrt(15)) / 37800
            W[GMS + 1, IW - 1] = (2665 - 14 * sqrt(15)) / 37800
            W[GMS + 2, IW - 1] = 2 * 15 / 567
            PTS[GMS + 3] = 0
        else:
            W[GMS, IW - 1] = (2 * (27 - N) / (N + 5) - L2 * (13 - N)) / (
                L1 ** 4 * (L1 - L2) * N4
            )
            W[GMS + 1, IW - 1] = (2 * (27 - N) / (N + 5) - L1 * (13 - N)) / (
                L2 ** 4 * (L2 - L1) * N4
            )
            W[GMS + 2, IW - 1] = (2 / (N + 5) - D2) / (N4 * (D1 - D2) * D1 ** 4)
            W[GMS + 3, IW - 1] = (2 / (N + 5) - D1) / (N4 * (D2 - D1) * D2 ** 4)
        IW = IW - 1
        G[0, 4] = 7 / (N + 7)
        G[1:NP, 4] = 1 / (N + 7)
        PTS[4] = NP
        G[0, 5] = 5 / (N + 7)
        G[1, 5] = 3 / (N + 7)
        G[2:NP, 5] = 1 / (N + 7)
        PTS[5] = NP * N
        G[0:3, 6] = 3 / (N + 7)
        if NP > 3:
            G[3:NP, 6] = 1 / (N + 7)
        PTS[6] = NP * N * (N - 1) / 6
        W[1, IW - 1] = (N + 3) ** 7 / (2 * 64 * N4 * (N + 5))
        W[2:4, IW - 1] = -((N + 5) ** 7) / (64 * N6)
        W[4:7, IW - 1] = (N + 7) ** 7 / (64 * N6 * (N + 7))
    if KEY == 4:
        IW = IW - 1
        SG = 1 / (23328 * N6)
        U5 = -(6 ** 3) * SG * (52212 - N * (6353 + N * (1934 - N * 27)))
        U6 = 6 ** 4 * SG * (7884 - N * (1541 - N * 9))
        U7 = -(6 ** 5) * SG * (8292 - N * (1139 - N * 3)) / (N + 7)
        P0 = -144 * (142528 + N * (23073 - N * 115))
        P1 = -12 * (6690556 + N * (2641189 + N * (245378 - N * 1495)))
        P2 = -16 * (6503401 + N * (4020794 + N * (787281 + N * (47323 - N * 385))))
        P3 = -(6386660 + N * (4411997 + N * (951821 + N * (61659 - N * 665)))) * (N + 7)
        A = P2 / (3 * P3)
        P = A * (P1 / P2 - A)
        Q = A * (2 * A * A - P1 / P3) + P0 / P3
        R = sqrt(-(P ** 3))
        TH = acos(-Q / (2 * R)) / 3
        R = 2 * R ** (1 / 3)
        TP = 2 * pi / 3
        A1 = -A + R * cos(TH)
        A2 = -A + R * cos(TH + 2 * TP)
        A3 = -A + R * cos(TH + TP)
        G[0, GMS + 4] = (1 - N * A1) / NP
        G[1:NP, GMS + 4] = (1 + A1) / NP
        PTS[GMS + 4] = NP
        G[0, GMS + 5] = (1 - N * A2) / NP
        G[1:NP, GMS + 5] = (1 + A2) / NP
        PTS[GMS + 5] = NP
        G[0, GMS + 6] = (1 - N * A3) / NP
        G[1:NP, GMS + 6] = (1 + A3) / NP
        PTS[GMS + 6] = NP
        W[GMS + 4, IW - 1] = (
            (U7 - (A2 + A3) * U6 + A2 * A3 * U5)
            / (A1 ** 2 - (A2 + A3) * A1 + A2 * A3)
            / A1 ** 5
        )
        W[GMS + 5, IW - 1] = (
            (U7 - (A1 + A3) * U6 + A1 * A3 * U5)
            / (A2 ** 2 - (A1 + A3) * A2 + A1 * A3)
            / A2 ** 5
        )
        W[GMS + 6, IW - 1] = (
            (U7 - (A2 + A1) * U6 + A2 * A1 * U5)
            / (A3 ** 2 - (A2 + A1) * A3 + A2 * A1)
            / A3 ** 5
        )
        G[0, GMS + 7] = 4 / (N + 7)
        G[1, GMS + 7] = 4 / (N + 7)
        G[2:NP, GMS + 7] = 1 / (N + 7)
        PTS[GMS + 7] = NP * N / 2
        W[GMS + 7, IW - 1] = 10 * (N + 7) ** 6 / (729 * N6)
        G[0, GMS + 8] = 11 / (N + 7) / 2
        G[1, GMS + 8] = 5 / (N + 7) / 2
        G[2:NP, GMS + 8] = 1 / (N + 7)
        PTS[GMS + 8] = NP * N
        W[GMS + 8, IW - 1] = 64 * (N + 7) ** 6 / (6561 * N6)
        W[3, IW - 1] = W[3, IW]
        W[6, IW - 1] = W[6, IW]
        IW = IW - 1
        G[0, 7] = 9 / (N + 9)
        G[1:NP, 7] = 1 / (N + 9)
        PTS[7] = NP
        G[0, 8] = 7 / (N + 9)
        G[1, 8] = 3 / (N + 9)
        G[2:NP, 8] = 1 / (N + 9)
        PTS[8] = NP * N
        G[0:2, 9] = 5 / (N + 9)
        G[2:NP, 9] = 1 / (N + 9)
        PTS[9] = NP * N / 2
        G[0, 10] = 5 / (N + 9)
        G[1:3, 10] = 3 / (N + 9)
        if NP > 3:
            G[3:NP, 10] = 1 / (N + 9)
        PTS[10] = NP * N * (N - 1) / 2
        W[1, IW - 1] = -((N + 3) ** 9) / (6 * 256 * N6)
        W[2:4, IW - 1] = (N + 5) ** 9 / (2 * 256 * N6 * (N + 7))
        W[4:7, IW - 1] = -((N + 7) ** 9) / (256 * N8)
        W[7:11, IW - 1] = (N + 9) ** 9 / (256 * N8 * (N + 9))
        if N > 2:
            G[0:4, 11] = 3 / (N + 9)
            if NP > 4:
                G[4:NP, 11] = 1 / (N + 9)
            PTS[11] = NP * N * (N - 1) * (N - 2) / 24
            W[11, IW - 1] = W[7, IW - 1]
    W[0, :] = 1 - PTS[1:WTS] @ W[1:WTS, :]
    NB = np.sum(PTS * (W[:, 0] * W[:, 0]))
    W[:, 1:RLS] = W[:, 1:RLS] - np.reshape(W[:, 0], (WTS, 1)) @ np.ones((1, RLS - 1))
    W[:, 1] = W[:, 1] * sqrt(NB / np.sum(PTS * W[:, 1] * W[:, 1]))
    for K in range(2, RLS):
        W[:, K] = (
            W[:, K]
            - W[:, 1:K] @ np.transpose(W[:, 1:K]) @ (PTS * W[:, K]) / NB
        )
        W[:, K] = W[:, K] * sqrt(NB / np.sum(PTS * W[:, K] * W[:, K]))
    return G, W, PTS


def _SMPRUL(ND, VRTS, VOL, NF, F, G, W, PTS):
    RTMN = 0.1
    SMALL = 1.0e-12
    ERRCOF = 8
    WTS, RLS = W.shape
    RULE = np.zeros((NF, RLS))
    for K in range(WTS):
        if PTS[K] > 0:
            RULE = RULE + VOL * np.reshape(
                _SMPSMS(ND, VRTS, NF, F, G[:, K]), (NF, 1)
            ) @ np.reshape(W[K, :], (1, RLS))
    BASVAL = np.zeros(NF)
    RGNERR = np.zeros(NF)
    for I in range(NF):
        BASVAL[I] = RULE[I, 0]
        NMBS = abs(BASVAL[I])
        RT = RTMN
        K = RLS - 1
        #NMCP = 1  ###
        while K >= 2:
            NMRL = max(abs(RULE[I, K]), abs(RULE[I, K - 1]))
            if (NMRL > SMALL * NMBS) & (K < RLS-1):
                RT = max(NMRL / NMCP, RT)
            RGNERR[I] = max(NMRL, RGNERR[I])
            NMCP = NMRL
            K = K - 2
        if (RT < 1) & (RLS > 2):
            RGNERR[I] = RT * NMCP
        RGNERR[I] = max(ERRCOF * RGNERR[I], SMALL * NMBS)
    return BASVAL, RGNERR


def _SMPSAD(ND, NF, F, MXFS, EA, ER, KEY, RCLS, SBS, VRTS, partitionInfo):
    NV = 0
    DFCOST = 1 + 2 * ND * (ND + 1)
    VL = np.zeros(NF)
    AE = np.zeros(NF)
    G, W, PTS = _SMPRMS(ND, KEY)
    FCT = factorial(ND)
    VLS = np.zeros(
        (NF, SBS)
    )  # VLS[i,j] = estimated integral of F[i] on simplex VRTS[,,j]
    AES = np.zeros(
        (NF, SBS)
    )  # AES[i,j] = estimated abs. err. in integral of F[i] on simplex VRTS[,,j]
    VOL = np.zeros(SBS)
    for K in range(SBS):
        VOL[K] = (
            abs(np.linalg.det(VRTS[K][:, 0:ND] - (VRTS[K][:, ND][:, np.newaxis]))) / FCT
        )
        BASVAL, RGNERR = _SMPRUL(ND, VRTS[K], VOL[K], NF, F, G, W, PTS)
        AES[:, K] = RGNERR
        VLS[:, K] = BASVAL
        VL = VL + VLS[:, K]
        AE = AE + AES[:, K]
        NV = NV + RCLS
    FL = int(any(np.greater(AE, np.max(np.concatenate(([EA], ER * np.abs(VL)))))))
    while (FL > 0) & (NV + DFCOST + 4 * RCLS <= MXFS):
        ID = np.argmax(np.apply_along_axis(np.max, 0, AES))
        VL = VL - VLS[:, ID]
        AE = AE - AES[:, ID]
        VRTS, NEW = _SMPDFS(ND, NF, F, ID, SBS, VRTS)
        VI = VOL[ID] / NEW
        VOL = np.concatenate((VOL, np.zeros(NEW - 1)))
        VLS = np.hstack((VLS, np.zeros((NF, NEW - 1))))
        AES = np.hstack((AES, np.zeros((NF, NEW - 1))))
        for K in np.concatenate(([ID], range(SBS, SBS + NEW - 1))):
            VOL[K] = VI
            BASVAL, RGNERR = _SMPRUL(ND, VRTS[K], VI, NF, F, G, W, PTS)
            VLS[:, K] = BASVAL
            AES[:, K] = RGNERR
            VL = VL + VLS[:, K]
            AE = AE + AES[:, K]
            NV = NV + RCLS
        NV = NV + DFCOST
        SBS = SBS + NEW - 1
        FL = int(any(np.greater(AE, np.max(np.concatenate(([EA], ER * np.abs(VL)))))))
    if SBS > 1:
        VL = np.apply_along_axis(np.sum, 1, VLS)
        AE = np.apply_along_axis(np.sum, 1, AES)
    if NF == 1:
        VL = VL[0]
        AE = AE[0]
    if partitionInfo:
        result = (VL, AE, NV, FL, VRTS, VLS, AES, VOL)
    else:
        result = (VL, AE, NV, FL)
    return result


def _adsimp(ND, VRTS, NF, F, MXFS, EA, ER, KEY, partitionInfo):
    if KEY == 0:
        KEY = 3
    SBS = len(VRTS)
    FL, RULCLS = _SMPCHC(ND, NF, MXFS, EA, ER, SBS, KEY)
    if FL != 0:
        return (np.repeat(None, NF), np.repeat(None, NF), 0, FL, None, None, None)
    return _SMPSAD(ND, NF, F, MXFS, EA, ER, KEY, RULCLS, SBS, VRTS, partitionInfo)


def _adsimp_message(rcode):
    if rcode == 0:
        msg = "OK"
    elif rcode == 1:
        msg = "error: maxEvals exceeded - too many function evaluations"
    elif rcode == 2:
        msg = "error: integRule < 0 or integRule > 4"
    elif rcode == 3:
        msg = "error: dimension n of the space < 2"
    elif rcode == 4:
        msg = "error: fDim = dimension of f < 1"
    elif rcode == 5:
        msg = "error: absError < 0 and tol < 0"
    else:
        msg = "error: unknown return code = " + str(rcode)
    return msg


def integrateOnSimplex(
    f, S, dim=1, maxEvals=10000, absError=0.0, tol=1.0e-5, rule=3, info=False
):
    """
    Integration on a simplex.
    
    Parameters
    ----------
    f : function
        The function to be integrated.
    S : array-like
        Simplex or simplices; a n-dimensional simplex is given as n+1 vectors of length n, the vertices.
    dim : integer
        The dimension of the values of `f`.
    maxEvals : integer
        Maximum number of calls to `f`.
    absError : number
        Desired absolute error.
    tol : number
        Desired relative error.
    rule : integer 
        Integer between 1 and 4 corresponding to the integration rule; a 2*rule+1 degree rule will be applied.
    info : Boolean
        Whether to return more info.
        
    Returns
    -------
    dictionary
        The value of the integral is in the field `"integral"` of the returned value.
        
    """
    S = np.asarray(S)
    if len(S.shape) == 2:
        S = np.array([S])
    nS = len(S)
    m, n = S[0].shape
    if m != n + 1:
        print("Invalid simplex")
        sys.exit(1)
    Simplices = np.empty((nS, n, m))
    for i in range(nS):
        if S[i].shape != (m, n):
            print("Invalid simplex found in `S`.")
            sys.exit(1)
        Simplices[i, :, :] = np.transpose(S[i])

    if dim == 1:
        fNew = lambda x : np.array([f(x)])
    else:
        fNew = lambda x : np.asarray(f(x))

    VL, AE, NV, FL = _adsimp(
        n, Simplices, dim, fNew, maxEvals, absError, tol, rule, info
    )
    result = {
        "integral": VL,
        "estAbsError": AE,
        "functionEvaluations": NV,
        "returnCode": FL,
        "message": _adsimp_message(FL),
    }
    return result


def _term(Q, monom):
    coef = Q.coeff_monomial(monom)
    powers = list(monom)
    j = sum(powers)
    if j == 0:
        return coef
    coef = coef * np.prod(list(map(factorial, powers)))
    n = len(monom)
    return coef / np.prod(list(range(n+1, n+j+1)))


def integratePolynomialOnSimplex(P, S):
    """
    Integration of a polynomial on a simplex.
    
    Parameters
    ----------
    P : polynomial
        The polynomial to be integrated.
    S : array-like
        Simplex; a n-dimensional simplex is given as n+1 vectors of length n, the vertices.
        
    Returns
    -------
    number
        The value of the integral of the polynomial over the simplex.
        
    """
    gens = P.gens
    n = len(gens)
    S = np.asarray(S)
    v = S[n,:]
    columns = []
    for i in range(n):
        columns.append(S[i,:] - v)    
    B = np.column_stack(tuple(columns))
    dico = {}
    for i in range(n):
        newvar = v[i]
        for j in range(n):
            newvar = newvar + B[i,j]*Poly(gens[j], gens, domain="RR")
        dico[gens[i]] = newvar.as_expr()
    Q = P.subs(dico, simultaneous=True).as_expr().as_poly(gens)
    monoms = Q.monoms()
    s = 0.0
    for monom in monoms:
        s = s + _term(Q, monom)
    return np.abs(np.linalg.det(B)) / factorial(n) * s
