import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from PIL import Image

favicon = Image.open("images/favicon.ico")

st.set_page_config(
    page_title="BD Hardening Soil",
    page_icon=favicon,
    layout="wide",
    initial_sidebar_state="expanded",
)

sb = st.sidebar
BD_im = Image.open("images/BDTunnelTools.png")
sb.image(BD_im)

with sb:
    st.caption(
        "You can find more on the theory of Hardening Soil in [the blog post here.](https://berkdemir.github.io/2021/05/11/Hardening-Soil-Model/)"
    )


st.header("Hardening Soil Model | Effect of Power m")

columns = st.columns([1, 1, 3])
with columns[0]:
    B = st.slider(
        "Width of Foundation (m)", min_value=1, max_value=100, step=1, value=20
    )
    L = st.slider(
        "Length of Foundation (m)", min_value=1, max_value=100, step=1, value=40
    )
    q = st.slider(
        "Pressure on Foundation (kPa)", min_value=0, max_value=500, step=10, value=100
    )
    UW = st.slider(
        "Unit Weight of Soil (kN/m3)", min_value=0, max_value=25, step=1, value=20
    )
    K0 = st.slider(
        "Coefficient of Lateral Earth Pressure",
        min_value=0.0,
        max_value=2.0,
        step=0.1,
        value=0.5,
    )

with columns[1]:
    m = st.slider("Power m", min_value=0.0, max_value=1.0, step=0.01, value=0.55)
    pref = st.slider(
        "Reference Pressure, pref (kPa)",
        min_value=10,
        max_value=500,
        step=10,
        value=100,
    )
    E50ref = st.slider("E50ref (MPa)", min_value=10, max_value=500, step=10, value=50)
    c = st.slider("Cohesion (kPa)", min_value=0, max_value=200, step=1, value=20)
    phi = st.slider(
        "Angle of Friction (deg)", min_value=0, max_value=45, step=1, value=30
    )
col2 = st.columns([2, 3])
with col2[-1]:
    max_depth = st.slider(
        "Maximum depth for figure (m)", min_value=10, max_value=100, value=20
    )


def boussinesq(L, B, z, q):
    m = L / B
    b = B / 2
    n = z / b

    I = (
        (
            m
            * n
            * (1 + m * m + 2 * n * n)
            / np.sqrt(1 + m * m + n * n)
            / (1 + n * n)
            / (m * m + n * n)
            + np.arcsin(m / (np.sqrt(m * m + n * n) * np.sqrt(1 + n * n)))
        )
        * 2
        / np.pi
    )
    return q * I


def Ecalc(E50ref, q, K0, m, pref, c, phi):
    return (
        E50ref
        * (
            (c * np.cos(np.radians(phi)) + q * K0 * np.sin(np.radians(phi)))
            / (c * np.cos(np.radians(phi)) + pref * np.sin(np.radians(phi)))
        )
        ** m
    )


depth = [i for i in range(0, max_depth + 1, 1)]
qred = list()
qred_noload = list()
E50_noload = list()
E50_load = list()
E50_load_m0 = list()
E50_load_m05 = list()
E50_load_m1 = list()
for i in depth:
    qred.append(boussinesq(L, B, i, q) + UW * i)
    E50_load.append(Ecalc(E50ref, qred[-1], K0, m, pref, c, phi))
    E50_load_m0.append(Ecalc(E50ref, qred[-1], K0, 0, pref, c, phi))
    E50_load_m05.append(Ecalc(E50ref, qred[-1], K0, 0.5, pref, c, phi))
    E50_load_m1.append(Ecalc(E50ref, qred[-1], K0, 1, pref, c, phi))
    print(E50_load_m1)


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6), sharey=True)
ax[0].plot(qred, depth, label="Pressure with Depth", color="red")
ax[0].set_xlabel("Pressure (kPa)")
ax[0].set_ylabel("Depth (m)")
ax[0].invert_yaxis()
ax[0].legend(fontsize=8)
ax[0].set_title("Pressure Difference with Depth (Boussinesq)")

ax[1].plot(E50_load, depth, label="m = " + str(round(m, 2)), color="red")
ax[1].plot(E50_load_m0, depth, label="m = 0", color="red", lw=1, ls="--")
ax[1].plot(E50_load_m05, depth, label="m = 0.5", color="green", lw=0.5, ls="--")
ax[1].plot(E50_load_m1, depth, label="m = 1", color="gray", lw=1, ls="--")

ax[1].set_xlabel("E50 (MPa)")
ax[1].legend(fontsize=8)
ax[1].set_title("E50 with Depth")
ax[1].fill_betweenx(depth, E50_load_m1, E50_load_m0, color="gray", alpha=0.1)

fig.suptitle("Effect of Load and Power m on E50", fontsize=16)
fig.tight_layout()


with columns[2]:
    st.pyplot(fig)