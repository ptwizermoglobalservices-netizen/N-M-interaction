import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

# ============================================
# USER INPUTS
# ============================================
st.title("RC Column N-M Interaction Diagram")

b = st.number_input("Section width b (mm)", value=300.0)
H = st.number_input("Section height H (mm)", value=500.0)
As = st.number_input("Tension reinforcement As (mm²)", value=2010.0)
Asp = st.number_input("Compression reinforcement As' (mm²)", value=804.0)

# ============================================
# DEFAULT VALUES
# ============================================
if st.button("Generate Interaction Diagram"):
    # Concrete
    fck = 30.0          # MPa
    alpha_cc = 0.85
    gamma_c = 1.5

    fcd = alpha_cc * fck / gamma_c      # MPa

    # Steel
    fyk = 450.0         # MPa
    gamma_s = 1.15

    fyd = fyk / gamma_s                 # MPa

    Es = 200000.0       # MPa

    ecu = 0.0035
    eyd = fyd / Es

    # stress block assumption
    xi = 0.809
    lam = 0.40

    # Reinforcement positions
    delta = 50.0        # Compression steel cover (mm)
    d = H - 50.0        # Effective depth (mm)

    # Section centroid
    Yg = H / 2.0

    # ============================================
    # STEEL STRESS FUNCTION
    # ============================================
    def steel_stress(eps):
        """
        Elastic-perfectly plastic steel model
        """
        if eps > eyd:
            return fyd
        elif eps < -eyd:
            return -fyd
        else:
            return Es * eps


    # ============================================
    # STORE RESULTS
    # ============================================
    results = []

    # ============================================
    # PURE TENSION
    # Xc -> -infinity
    # ============================================
    Nu = -fyd * (As + Asp)
    Mu = 0.0

    results.append({
        'Xc': -np.inf,
        'sigma_s': -fyd,
        "sigma_s'": -fyd,
        'C': 0.0,
        "C'": -fyd * Asp,
        'T': fyd * As,
        'Nu (KN)': Nu/1000,
        'Mu (kN.m)': Mu
    })

    # ============================================
    # 0 <= Xc <= H
    # ============================================
    Xc_values = np.linspace(1, H, 200)

    for Xc in Xc_values:

        # Steel strains
        eps_sp = ((Xc - delta) / Xc) * ecu
        eps_s = ((d - Xc) / Xc) * ecu

        # Steel stresses
        sigma_sp = steel_stress(eps_sp)
        sigma_s = steel_stress(eps_s)

        # Forces
        C = xi * fcd * b * Xc
        Cp = sigma_sp * Asp
        T = sigma_s * As

        # Axial force
        Nu = C + Cp - T

        # Moment about centroid
        Mu = (
            C * (Yg - lam * Xc)
            + Cp * (Yg - delta)
            + T * (d - Yg)
        )

        Mu = Mu / 1e6      # N.mm → kN.m

        results.append({
            'Xc': Xc,
            'sigma_s': sigma_s,
            "sigma_s'": sigma_sp,
            'C': C,
            "C'": Cp,
            'T': T,
            'Nu (KN)': Nu/1000,
            'Mu (kN.m)': Mu
        })

    # ============================================
    # PURE COMPRESSION
    # Xc -> +infinity
    # Entire concrete uniformly stressed
    # ============================================
    Nu = b * H * fcd + fyd * (As + Asp)
    Mu = 0.0

    results.append({
        'Xc': np.inf,
        'sigma_s': fyd,
        "sigma_s'": fyd,
        'C': b * H * fcd,
        "C'": fyd * Asp,
        'T': -fyd * As,
        'Nu (KN)': Nu/1000,
        'Mu (kN.m)': Mu
    })

    # ============================================
    # CREATE TABLE
    # ============================================
    df = pd.DataFrame(results)

    print("\nInteraction Table:")
    st.dataframe(df)

    # ============================================
    # SAVE RESULTS
    # ============================================
    df.to_csv("NM_interaction_results.csv", index=False)

    print("\nResults saved as:")
    print("NM_interaction_results.csv")

    # ============================================
    # PLOT N-M CURVE
    # ============================================
    fig =plt.figure(figsize=(8,6))

    plt.plot(
        df['Mu (kN.m)'],
        df['Nu (KN)'],
        '-o',
        markersize=3
    )

    plt.xlabel('Moment Mu (kN.m)')
    plt.ylabel('Axial Force Nu (N)')
    plt.title('N-M Interaction Diagram')

    plt.grid(True)

    st.pyplot(fig)
