"""
DSAS-style Shoreline Change Analysis — Cahuita / Gandoca-Manzanillo
====================================================================
Implementa fielmente la metodología del USGS DSAS v5.1 (Himmelstoss et al., 2021).


Tasas:
  EPR  – End Point Rate            [m/año]  + incertidumbre 95%
  LRR  – Linear Regression Rate    [m/año]  + CI95, SE, R², p-valor
  NSM  – Net Shoreline Movement    [metros]
  SCE  – Shoreline Change Envelope [metros]
  WLR  – Weighted Linear Regression [m/año] (ponderada por 1/σ²)

Capas de salida:
  1. DSAS_EPR           — transectas coloreadas por EPR (9 bins rojo-azul)
  2. DSAS_LRR           — ídem por LRR
  3. DSAS_NSM           — ídem por NSM
  4. DSAS_SCE           — ídem por SCE
  5. DSAS_WLR           — ídem por WLR
  6. DSAS_Forecast_Zona — polígono incertidumbre 95% a N años
  7. DSAS_Forecast_Lin  — línea central del forecast (LRR)

USO: QGIS → Complementos → Consola Python → Mostrar Editor → ▶
"""

# ─────────────────────────────────────────────────────────────────────────────
#  DEPENDENCIAS
# ─────────────────────────────────────────────────────────────────────────────
from qgis.core import (
    QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry,
    QgsPointXY, QgsField, QgsWkbTypes,
    QgsGraduatedSymbolRenderer, QgsRendererRange,
    QgsVectorFileWriter, QgsLineSymbol, QgsFillSymbol,
)
from qgis.PyQt.QtCore import QVariant
import math, datetime, os
from collections import defaultdict


# ═════════════════════════════════════════════════════════════════════════════
#  CONFIGURACIÓN
# ═════════════════════════════════════════════════════════════════════════════

NOMBRE_CAPA_SHORELINES  = "shorelines_shape2"
CAMPO_FECHA             = "FECHA"          # campo de fecha en la capa (texto)
NOMBRE_BASELINE         = "BASELINE onshore"

ESPACIADO_TRANSECTAS    = 50    # metros entre transectas
LARGO_TRANSECTA         = 250   # metros de búsqueda desde el baseline (cada lado)

# ── Suavizado de la dirección perpendicular (igual que DSAS "smoothing distance")
# DSAS recomienda un valor mayor que el ancho de las curvas de la costa.
# Para Cahuita/Gandoca: 200–500 m suele funcionar bien.
SMOOTHING_DIST          = 300   # metros

# ── Orientación tierra/mar (CRÍTICO — igual que "Identify land location" en DSAS)
# "auto" : el script auto-detecta usando las primeras 20 transectas
# "izq"  : tierra a la IZQUIERDA del baseline según dirección de digitalización
# "der"  : tierra a la DERECHA del baseline según dirección de digitalización
TIERRA_ORIENTACION      = "auto"

# ── Intersección a usar cuando la transecta corta la misma shoreline >1 vez
# "seaward"  → la más alejada del baseline hacia el mar   ← recomendado (como DSAS)
# "landward" → la más cercana al baseline
INTERSECCION            = "seaward"

# ── Incertidumbre posicional por imagen (metros)
# Sentinel-2 ~5 m  |  Landsat ~15 m
# DSAS usa este valor en el campo DSAS_uncy de cada shoreline
INCERTIDUMBRE_DEFAULT   = 5.0
INCERTIDUMBRE_FECHAS    = {
    # "01/06/2002": 15.0,
    # "01/06/2024": 5.0,
}

# ── Forecast
AÑOS_FORECAST           = 10   # años a proyectar

# ── Salida
CARPETA_SALIDA          = ""   # dejar vacío = solo capas temporales en QGIS, sin guardar archivos


# ═════════════════════════════════════════════════════════════════════════════
#  MATEMÁTICAS  (idénticas a DSAS v5.1)
# ═════════════════════════════════════════════════════════════════════════════

_T95 = {1:12.706,2:4.303,3:3.182,4:2.776,5:2.571,6:2.447,7:2.365,
        8:2.306,9:2.262,10:2.228,15:2.131,20:2.086,30:2.042,60:2.000,120:1.980}

def t_critico(df):
    if df <= 0: return float('inf')
    for k in sorted(_T95):
        if df <= k: return _T95[k]
    return 1.960

def p_aprox(t_stat, df):
    if df < 3 or t_stat is None: return None
    t2 = t_stat**2
    z = abs(t_stat)*(1 - 1/(4*df)) / math.sqrt(1 + t2/(2*df))
    x = z / math.sqrt(2)
    a = [0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429]
    p = 0.3275911
    tv = 1/(1+p*x)
    y = 1-(((((a[4]*tv+a[3])*tv)+a[2])*tv+a[1])*tv+a[0])*tv*math.exp(-x*x)
    return round(2*(1-y), 4)


def regresion_ols(xs, ys, ws=None):
    """
    OLS o WLS (si ws provisto).  Devuelve dict completo estilo DSAS.
    xs = años decimales, ys = distancias firmadas, ws = pesos 1/σ²
    """
    n = len(xs)
    if n < 2: return None
    if ws is None: ws = [1.0]*n
    W  = sum(ws)
    Wx = sum(w*x   for w,x in zip(ws,xs))
    Wy = sum(w*y   for w,y in zip(ws,ys))
    Wx2= sum(w*x*x for w,x in zip(ws,xs))
    Wxy= sum(w*x*y for w,x,y in zip(ws,xs,ys))
    Sxx= Wx2 - Wx**2/W
    Sxy= Wxy - Wx*Wy/W
    if Sxx == 0: return None
    slope     = Sxy/Sxx
    intercept = (Wy - slope*Wx)/W
    x_mean    = Wx/W
    y_pred    = [slope*x+intercept for x in xs]
    ss_res    = sum(w*(y-yp)**2 for w,y,yp in zip(ws,ys,y_pred))
    ss_tot    = sum(w*(y-Wy/W)**2 for w,y in zip(ws,ys))
    r2        = 1-ss_res/ss_tot if ss_tot>0 else 1.0
    if n > 2:
        mse = ss_res/(n-2)
        se  = math.sqrt(mse/Sxx) if Sxx>0 else 0.0
        tc  = t_critico(n-2)
        lci, uci = slope-tc*se, slope+tc*se
        t_st = slope/se if se>0 else None
        pval = p_aprox(t_st, n-2)
    else:
        mse=se=lci=uci=t_st=pval=None
    return dict(slope=slope, intercept=intercept, x_mean=x_mean, Sxx=Sxx,
                r2=r2, se=se, lci=lci, uci=uci, n=n, mse=mse, pval=pval)


def prediccion_pi(res, x_new):
    """Intervalo de predicción al 95% en x_new. Devuelve (central, lci, uci)."""
    pos = res['slope']*x_new + res['intercept']
    if res['n'] < 3 or res['mse'] is None: return pos, None, None
    se_pred = math.sqrt(res['mse']*(1+1/res['n']+(x_new-res['x_mean'])**2/res['Sxx']))
    tc      = t_critico(res['n']-2)
    m       = tc*se_pred
    return pos, pos-m, pos+m


def parsear_fecha(valor):
    """DD/MM/YYYY, YYYY-MM-DD, etc. → año decimal float."""
    s = str(valor).strip()
    for fmt in ("%d/%m/%Y","%Y-%m-%d","%Y/%m/%d","%d-%m-%Y"):
        try:
            dt   = datetime.datetime.strptime(s,fmt)
            base = datetime.datetime(dt.year,1,1)
            sig  = datetime.datetime(dt.year+1,1,1)
            return dt.year + (dt-base).days/(sig-base).days
        except ValueError: continue
    try: return float(s[:4])
    except: return None


def fusionar_geom(capa):
    g = QgsGeometry()
    for f in capa.getFeatures():
        gg = f.geometry()
        if gg and not gg.isEmpty():
            g = gg if g.isEmpty() else g.combine(gg)
    return g


def primera_linea(geom):
    if geom.isMultipart():
        p = geom.asGeometryCollection()
        return p[0] if p else geom
    return geom


def todos_los_puntos(geom):
    """Extrae TODOS los puntos de cualquier geometría de intersección."""
    pts = []
    if geom.isEmpty(): return pts
    t = QgsWkbTypes.flatType(geom.wkbType())
    if t == QgsWkbTypes.Point:
        pts.append(geom.asPoint())
    elif t == QgsWkbTypes.MultiPoint:
        pts.extend(geom.asMultiPoint())
    elif t in (QgsWkbTypes.GeometryCollection, QgsWkbTypes.LineString,
               QgsWkbTypes.MultiLineString):
        for g in geom.asGeometryCollection():
            pts.extend(todos_los_puntos(g))
    return pts



# ═════════════════════════════════════════════════════════════════════════════
#  RENDERERS ESTILO DSAS (9 bins: 4 rojo · 1 gris · 4 azul)
# ═════════════════════════════════════════════════════════════════════════════

def renderer_dsas_9bins(campo, valores):
    """
    Reproduce exactamente la rampa de DSAS:
    4 tonos rojos (negativo) · gris neutro (≈0) · 4 tonos azules (positivo).
    Los cortes se calculan al percentil 85 de los datos, como hace DSAS.
    """
    vals = sorted(v for v in valores if v is not None)
    if not vals:
        return QgsGraduatedSymbolRenderer(campo, [])

    neg_vals = sorted(abs(v) for v in vals if v < 0)
    pos_vals = sorted(v for v in vals if v > 0)

    # Percentil 85 para cada lado
    def p85(lst):
        if not lst: return 1.0
        return lst[int(len(lst)*0.85)]

    p85_neg = p85(neg_vals) if neg_vals else 1.0
    p85_pos = p85(pos_vals) if pos_vals else 1.0

    # 4 bins por lado + 1 central
    step_neg = p85_neg / 4
    step_pos = p85_pos / 4
    neutral  = min(step_neg, step_pos) * 0.5   # zona central ±neutral

    colores_neg = ["#67001f","#b2182b","#d6604d","#f4a582"]  # oscuro→claro
    colores_pos = ["#92c5de","#4393c3","#2166ac","#053061"]  # claro→oscuro

    rangos = []

    # Negativos (erosión) — de más extremo a zona central
    # El primer límite usa min(vals)-1 para capturar TODOS los valores más extremos que p85
    limites_neg = [min(vals) - 1, -p85_neg,
                   -p85_neg*2/3, -p85_neg/3, -neutral]
    for i in range(4):
        lo = limites_neg[i]
        hi = limites_neg[i+1]
        if hi < min(vals): continue
        sym = QgsLineSymbol.createSimple({"color": colores_neg[i], "width":"0.9"})
        rangos.append(QgsRendererRange(lo, hi, sym,
                      f"{lo:.2f} a {hi:.2f} m"))

    # Central (gris)
    sym_neu = QgsLineSymbol.createSimple({"color":"#888888","width":"0.7"})
    rangos.append(QgsRendererRange(-neutral, neutral, sym_neu,
                  f"-{neutral:.2f} a +{neutral:.2f} m (estable)"))

    # Positivos (acreción)
    # El último límite usa max(vals)+1 para capturar TODOS los valores más extremos que p85
    limites_pos = [neutral, p85_pos/3, p85_pos*2/3,
                   p85_pos, max(vals) + 1]
    for i in range(4):
        lo = limites_pos[i]
        hi = limites_pos[i+1]
        if lo > max(vals): continue
        sym = QgsLineSymbol.createSimple({"color": colores_pos[i], "width":"0.9"})
        rangos.append(QgsRendererRange(lo, hi, sym,
                      f"{lo:.2f} a {hi:.2f} m"))

    r = QgsGraduatedSymbolRenderer(campo, rangos)
    return r


# ═════════════════════════════════════════════════════════════════════════════
#  CARGA DE CAPAS
# ═════════════════════════════════════════════════════════════════════════════

print("="*65)
print("  DSAS v5.1-style  —  Shoreline Change Analysis")
print("="*65)

def get_capa(nombre):
    capas = QgsProject.instance().mapLayersByName(nombre)
    if not capas:
        raise ValueError(f"❌ Capa '{nombre}' no encontrada en el panel.")
    return capas[0]

baseline_capa   = get_capa(NOMBRE_BASELINE)
shorelines_capa = get_capa(NOMBRE_CAPA_SHORELINES)
baseline_geom   = primera_linea(fusionar_geom(baseline_capa))
largo_total     = baseline_geom.length()
crs_str         = baseline_capa.crs().authid()

# ── Indexar shorelines ─────────────────────────────────────────────────────
campos_disp = [f.name() for f in shorelines_capa.fields()]
if CAMPO_FECHA not in campos_disp:
    raise ValueError(f"❌ Campo '{CAMPO_FECHA}' no existe. Disponibles: {campos_disp}")

print(f"\nIndexando shorelines (campo '{CAMPO_FECHA}')...")
geoms_dict = defaultdict(list)
raw_fecha  = {}
uncy_dict  = {}

for feat in shorelines_capa.getFeatures():
    val = feat[CAMPO_FECHA]
    ad  = parsear_fecha(val)
    if ad is None:
        print(f"  ⚠ Fecha no parseable: '{val}' — ignorada")
        continue
    geoms_dict[ad].append(feat.geometry())
    raw_fecha[ad] = str(val).strip()
    uncy_dict[ad] = INCERTIDUMBRE_FECHAS.get(str(val).strip(), INCERTIDUMBRE_DEFAULT)

geoms_finales = {}
for ad, geoms in geoms_dict.items():
    union = QgsGeometry()
    for g in geoms:
        if g and not g.isEmpty():
            union = g if union.isEmpty() else union.combine(g)
    if not union.isEmpty():
        geoms_finales[ad] = union

fechas_ord = sorted(geoms_finales.keys())
n_fechas   = len(fechas_ord)
if n_fechas < 2:
    raise ValueError(f"❌ Se necesitan ≥2 fechas; encontradas: {n_fechas}")

ad_min, ad_max = fechas_ord[0], fechas_ord[-1]
span           = ad_max - ad_min
ad_forecast    = ad_max + AÑOS_FORECAST

print(f"  {n_fechas} shorelines cargadas:")
for ad in fechas_ord:
    print(f"    {raw_fecha[ad]}  (±{uncy_dict[ad]:.0f} m)")
print(f"\n  Período: {raw_fecha[ad_min]} → {raw_fecha[ad_max]}  ({span:.1f} años)")
print(f"  Forecast: año {int(ad_forecast)} (+{AÑOS_FORECAST} a)")


# ═════════════════════════════════════════════════════════════════════════════
#  AUTO-DETECCIÓN DE ORIENTACIÓN TIERRA/MAR
# ═════════════════════════════════════════════════════════════════════════════

def calcular_perp(dist, flip=False):
    """Perpendicular al baseline en la posición 'dist', con suavizado DSAS."""
    d1 = max(0.0, dist - SMOOTHING_DIST/2)
    d2 = min(largo_total, dist + SMOOTHING_DIST/2)
    p1 = baseline_geom.interpolate(d1).asPoint()
    p2 = baseline_geom.interpolate(d2).asPoint()
    dx, dy = p2.x()-p1.x(), p2.y()-p1.y()
    lon = math.hypot(dx, dy)
    if lon == 0: return None, None
    # Rotación 90° antihoraria = perpendicular izquierda
    px, py = -dy/lon, dx/lon
    if flip: px, py = -px, -py
    return px, py

if TIERRA_ORIENTACION == "auto":
    print("\nAuto-detección de orientación tierra/mar...")
    n_test      = min(20, int(largo_total/ESPACIADO_TRANSECTAS))
    positivos_n = 0
    negativos_n = 0
    for i in range(n_test):
        dist  = i * ESPACIADO_TRANSECTAS
        px,py = calcular_perp(dist, flip=False)
        if px is None: continue
        punto_base = baseline_geom.interpolate(dist).asPoint()
        for ad, gsh in list(geoms_finales.items())[:3]:
            t_i = QgsPointXY(punto_base.x()-px*LARGO_TRANSECTA,
                             punto_base.y()-py*LARGO_TRANSECTA)
            t_f = QgsPointXY(punto_base.x()+px*LARGO_TRANSECTA,
                             punto_base.y()+py*LARGO_TRANSECTA)
            transecta = QgsGeometry.fromPolylineXY([t_i, t_f])
            inter = transecta.intersection(gsh)
            for pt in todos_los_puntos(inter):
                d = (pt.x()-punto_base.x())*px + (pt.y()-punto_base.y())*py
                if d > 0.5: positivos_n += 1
                elif d < -0.5: negativos_n += 1
    FLIP = (negativos_n > positivos_n)
    dire = "der→izq (flip)" if FLIP else "izq→der (normal)"
    print(f"  Positivos: {positivos_n} | Negativos: {negativos_n} → {dire}")
elif TIERRA_ORIENTACION == "der":
    FLIP = True
    print("Orientación: tierra a la DERECHA del baseline (flip)")
else:
    FLIP = False
    print("Orientación: tierra a la IZQUIERDA del baseline (normal)")


# ═════════════════════════════════════════════════════════════════════════════
#  DEFINICIÓN DE CAPAS DE SALIDA
# ═════════════════════════════════════════════════════════════════════════════

def nueva_capa_linea(nombre, campos):
    c = QgsVectorLayer(f"LineString?crs={crs_str}", nombre, "memory")
    pr = c.dataProvider(); pr.addAttributes(campos); c.updateFields()
    return c, pr

campos_epr = [QgsField("TransectID",QVariant.Int), QgsField("EPR",QVariant.Double),
              QgsField("EPRunc",QVariant.Double), QgsField("EPR_EY",QVariant.Double),
              QgsField("EPR_OY",QVariant.Double), QgsField("EPR_Yrs",QVariant.Double)]

campos_lrr = [QgsField("TransectID",QVariant.Int), QgsField("LRR",QVariant.Double),
              QgsField("LCI95",QVariant.Double), QgsField("UCI95",QVariant.Double),
              QgsField("LSE",QVariant.Double), QgsField("LR2",QVariant.Double),
              QgsField("Lp_val",QVariant.Double), QgsField("Ldn",QVariant.Int)]

campos_nsm = [QgsField("TransectID",QVariant.Int), QgsField("NSM",QVariant.Double),
              QgsField("NSM_ODate",QVariant.String), QgsField("NSM_EDate",QVariant.String),
              QgsField("NSM_Yrs",QVariant.Double)]

campos_sce = [QgsField("TransectID",QVariant.Int), QgsField("SCE",QVariant.Double),
              QgsField("SCE_DDate",QVariant.String), QgsField("SCE_ADate",QVariant.String)]

campos_wlr = [QgsField("TransectID",QVariant.Int), QgsField("WLR",QVariant.Double),
              QgsField("WCI95",QVariant.Double), QgsField("WUI95",QVariant.Double),
              QgsField("WSE",QVariant.Double), QgsField("WR2",QVariant.Double),
              QgsField("Wdn",QVariant.Int)]

campos_fore = [QgsField("TransectID",QVariant.Int), QgsField("FORE_Yr",QVariant.Int),
               QgsField("FORE_Pos",QVariant.Double), QgsField("FORE_LCI",QVariant.Double),
               QgsField("FORE_UCI",QVariant.Double), QgsField("FORE_LRR",QVariant.Double)]

capa_epr, pr_epr = nueva_capa_linea("DSAS_EPR", campos_epr)
capa_lrr, pr_lrr = nueva_capa_linea("DSAS_LRR", campos_lrr)
capa_nsm, pr_nsm = nueva_capa_linea("DSAS_NSM", campos_nsm)
capa_sce, pr_sce = nueva_capa_linea("DSAS_SCE", campos_sce)
capa_wlr, pr_wlr = nueva_capa_linea("DSAS_WLR", campos_wlr)
capa_fore_lin = QgsVectorLayer(f"LineString?crs={crs_str}", "DSAS_Forecast_Lin",  "memory")
capa_fore_zon = QgsVectorLayer(f"Polygon?crs={crs_str}",   "DSAS_Forecast_Zona", "memory")
pr_fore_lin = capa_fore_lin.dataProvider()
pr_fore_zon = capa_fore_zon.dataProvider()
pr_fore_lin.addAttributes(campos_fore);  capa_fore_lin.updateFields()
pr_fore_zon.addAttributes([QgsField("FORE_Yr",QVariant.Int), QgsField("Conf",QVariant.String),
                            QgsField("n_trans",QVariant.Int)]); capa_fore_zon.updateFields()


# ═════════════════════════════════════════════════════════════════════════════
#  BUCLE PRINCIPAL: TRANSECTAS
# ═════════════════════════════════════════════════════════════════════════════

n_est = int(largo_total / ESPACIADO_TRANSECTAS)
print(f"\nGenerando {n_est} transectas cada {ESPACIADO_TRANSECTAS} m (±{LARGO_TRANSECTA} m, smooth {SMOOTHING_DIST} m)...")

feats_epr,feats_lrr,feats_nsm,feats_sce,feats_wlr = [],[],[],[],[]
fore_data  = []   # datos para post-proceso del forecast ANTES de armar el polígono
sin_datos  = 0
mult_inter = 0    # transectas donde se descartó alguna intersección múltiple

for i in range(n_est+1):
    dist = i * ESPACIADO_TRANSECTAS
    if dist > largo_total: break

    # ── Geometría de la transecta (perpendicular suavizada) ────────────────
    punto_base = baseline_geom.interpolate(dist).asPoint()
    px, py = calcular_perp(dist, flip=FLIP)
    if px is None: continue

    t_ini = QgsPointXY(punto_base.x()-px*LARGO_TRANSECTA,
                       punto_base.y()-py*LARGO_TRANSECTA)
    t_fin = QgsPointXY(punto_base.x()+px*LARGO_TRANSECTA,
                       punto_base.y()+py*LARGO_TRANSECTA)
    transecta = QgsGeometry.fromPolylineXY([t_ini, t_fin])

    def dfirmada(pt):
        return (pt.x()-punto_base.x())*px + (pt.y()-punto_base.y())*py

    # ── Intersecciones con cada shoreline ─────────────────────────────────
    # DSAS: selecciona intersección más seaward (o landward) cuando hay varias
    posiciones = {}
    for ad, gsh in geoms_finales.items():
        inter = transecta.intersection(gsh)
        pts   = todos_los_puntos(inter)
        if not pts: continue

        dists_pts = [(dfirmada(pt), pt) for pt in pts]
        # Filtrar: solo lado seaward (positivo), igual que baseline onshore en DSAS
        pos_positivos = [(d, pt) for d, pt in dists_pts if d > 0.5]

        if pos_positivos:
            if len(pos_positivos) > 1: mult_inter += 1
            if INTERSECCION == "seaward":
                d_sel, _ = max(pos_positivos, key=lambda x: x[0])
            else:
                d_sel, _ = min(pos_positivos, key=lambda x: x[0])
            posiciones[ad] = d_sel
        elif dists_pts:
            # Ninguna positiva → baseline podría ser offshore o error
            # Tomar la menos negativa
            d_sel, _ = max(dists_pts, key=lambda x: x[0])
            posiciones[ad] = d_sel

    n_pos = len(posiciones)
    if n_pos < 2:
        sin_datos += 1
        continue

    fechas_con = sorted(posiciones.keys())
    f_ant, f_nue = fechas_con[0], fechas_con[-1]
    span_t = f_nue - f_ant
    if span_t <= 0:
        sin_datos += 1
        continue

    vals = list(posiciones.values())
    xs   = list(posiciones.keys())
    ys   = [posiciones[a] for a in xs]

    # ── EPR  (DSAS: oldest→newest / years) ────────────────────────────────
    nsm_val = posiciones[f_nue] - posiciones[f_ant]
    epr_val = nsm_val / span_t
    # Incertidumbre EPR: √(σ_ant² + σ_nue²) / Δt  (Himmelstoss et al. 2021, §7.6.1)
    ua_ant  = uncy_dict[f_ant]
    ua_nue  = uncy_dict[f_nue]
    epr_unc = math.sqrt(ua_ant**2 + ua_nue**2) / span_t

    feat = QgsFeature(); feat.setGeometry(transecta)
    feat.setAttributes([i, round(epr_val,3), round(epr_unc,3),
                        round(f_nue,2), round(f_ant,2), round(span_t,2)])
    feats_epr.append(feat)

    # ── LRR  (OLS, §7.4) ─────────────────────────────────────────────────
    res_lrr = regresion_ols(xs, ys)
    lrr_val = res_lrr['slope'] if res_lrr else epr_val

    if res_lrr:
        feat = QgsFeature(); feat.setGeometry(transecta)
        feat.setAttributes([i, round(lrr_val,3),
                             round(res_lrr['lci'],3) if res_lrr['lci'] is not None else None,
                             round(res_lrr['uci'],3) if res_lrr['uci'] is not None else None,
                             round(res_lrr['se'],4)  if res_lrr['se']  is not None else None,
                             round(res_lrr['r2'],4),
                             round(res_lrr['pval'],4) if res_lrr['pval'] is not None else None,
                             res_lrr['n']])
        feats_lrr.append(feat)

    # ── WLR  (OLS ponderada 1/σ², §7.5) ──────────────────────────────────
    ws_raw = [1.0/(uncy_dict[a]**2) for a in xs]
    res_wlr = regresion_ols(xs, ys, ws=ws_raw)
    wlr_val = res_wlr['slope'] if res_wlr else lrr_val

    if res_wlr:
        feat = QgsFeature(); feat.setGeometry(transecta)
        feat.setAttributes([i, round(wlr_val,3),
                             round(res_wlr['lci'],3) if res_wlr['lci'] is not None else None,
                             round(res_wlr['uci'],3) if res_wlr['uci'] is not None else None,
                             round(res_wlr['se'],4)  if res_wlr['se']  is not None else None,
                             round(res_wlr['r2'],4),
                             res_wlr['n']])
        feats_wlr.append(feat)

    # ── NSM ───────────────────────────────────────────────────────────────
    feat = QgsFeature(); feat.setGeometry(transecta)
    feat.setAttributes([i, round(nsm_val,2), raw_fecha[f_ant], raw_fecha[f_nue],
                        round(span_t,2)])
    feats_nsm.append(feat)

    # ── SCE ───────────────────────────────────────────────────────────────
    sce_val  = max(vals) - min(vals)
    ad_max_v = max(posiciones, key=posiciones.get)
    ad_min_v = min(posiciones, key=posiciones.get)
    feat = QgsFeature(); feat.setGeometry(transecta)
    feat.setAttributes([i, round(sce_val,2), raw_fecha[ad_max_v], raw_fecha[ad_min_v]])
    feats_sce.append(feat)

    # ── FORECAST (LRR) — recopilar datos; el polígono se arma DESPUÉS del bucle ──
    # Solo se incluyen transectas con n ≥ 3 (mínimo para calcular PI)
    if res_lrr and res_lrr['n'] >= 3 and res_lrr['mse'] is not None:
        f_pos, f_lci, f_uci = prediccion_pi(res_lrr, ad_forecast)
        if f_lci is not None and f_uci is not None:
            fore_data.append({
                'pb': punto_base, 'px': px, 'py': py,
                'central': f_pos, 'uci': f_uci, 'lci': f_lci,
                'pi_width': f_uci - f_lci,
                'lrr': lrr_val, 'n': res_lrr['n'], 'tid': i,
            })


# ── Agregar features ────────────────────────────────────────────────────────
pr_epr.addFeatures(feats_epr); capa_epr.updateExtents()
pr_lrr.addFeatures(feats_lrr); capa_lrr.updateExtents()
pr_nsm.addFeatures(feats_nsm); capa_nsm.updateExtents()
pr_sce.addFeatures(feats_sce); capa_sce.updateExtents()
pr_wlr.addFeatures(feats_wlr); capa_wlr.updateExtents()

total = len(feats_epr)
print(f"  ✓ {total} transectas válidas  |  {sin_datos} sin intersección  |  {mult_inter} con intersección múltiple resuelta")


# ═════════════════════════════════════════════════════════════════════════════
#  POST-PROCESO DEL FORECAST: filtrar, acotar y suavizar
# ═════════════════════════════════════════════════════════════════════════════

print(f"\nPost-proceso forecast: {len(fore_data)} transectas con PI calculable...")

if len(fore_data) >= 3:

    # ── 1. Acotar PI extremos (percentil 95 del ancho) ────────────────────
    # Estadísticamente: un PI muy ancho indica MSE extremo (outlier local).
    # Se acota al p95 × 1.5 (regla similar a los bigotes de Tukey).
    # El valor central del forecast NO se cambia, solo el ancho del intervalo.
    anchos = sorted(d['pi_width'] for d in fore_data)
    n_fd   = len(anchos)
    p95_w  = anchos[min(int(n_fd*0.95), n_fd-1)]
    ancho_max = p95_w * 1.5
    n_acotados = 0
    for d in fore_data:
        if d['pi_width'] > ancho_max:
            mitad     = ancho_max / 2
            d['uci']  = d['central'] + mitad
            d['lci']  = d['central'] - mitad
            d['pi_width'] = ancho_max
            n_acotados += 1
    print(f"  PI p95={p95_w:.1f} m  |  máx permitido={ancho_max:.1f} m  |  {n_acotados} acotados")

    # ── 2. Suavizado geográfico de los límites (ventana deslizante) ───────
    # Igual a lo que hace DSAS internamente para producir un polígono limpio.
    # El suavizado actúa sobre las COORDENADAS de los puntos del polígono,
    # no sobre los valores estadísticos de la tabla de atributos.
    # Ventana W=7 → cada punto promedia ±3 vecinos.
    SMOOTH_W = 7

    def suavizar_lista(vals):
        n = len(vals); out = []
        hw = SMOOTH_W // 2
        for j in range(n):
            lo, hi = max(0, j-hw), min(n, j+hw+1)
            out.append(sum(vals[lo:hi]) / (hi-lo))
        return out

    # Extraer coordenadas de upper y lower en el espacio geográfico
    def pts_coords(key):
        return [(d['pb'].x() + d['px']*d[key],
                 d['pb'].y() + d['py']*d[key]) for d in fore_data]

    cent_coords  = pts_coords('central')
    upper_coords = pts_coords('uci')
    lower_coords = pts_coords('lci')

    # Suavizar X e Y por separado
    def smooth_xy(coords):
        xs = suavizar_lista([c[0] for c in coords])
        ys = suavizar_lista([c[1] for c in coords])
        return [QgsPointXY(x, y) for x, y in zip(xs, ys)]

    pts_central_sm = smooth_xy(cent_coords)
    pts_upper_sm   = smooth_xy(upper_coords)
    pts_lower_sm   = smooth_xy(lower_coords)
    print(f"  Suavizado W={SMOOTH_W}: ✓")

    # ── 3. Línea central del forecast ─────────────────────────────────────
    lrr_vals_all = [f["LRR"] for f in capa_lrr.getFeatures() if f["LRR"] is not None]
    lrr_prom     = sum(lrr_vals_all)/len(lrr_vals_all) if lrr_vals_all else 0

    feat_lin = QgsFeature()
    feat_lin.setGeometry(QgsGeometry.fromPolylineXY(pts_central_sm))
    feat_lin.setAttributes([0, int(ad_forecast), None, None, None, round(lrr_prom,3)])
    pr_fore_lin.addFeatures([feat_lin]); capa_fore_lin.updateExtents()

    # ── 4. Polígono de incertidumbre suavizado ────────────────────────────
    # Borde superior (seaward) → borde inferior (landward) invertido → cierre
    anillo = pts_upper_sm + list(reversed(pts_lower_sm)) + [pts_upper_sm[0]]
    geom_p = QgsGeometry.fromPolygonXY([anillo])
    if geom_p and not geom_p.isEmpty():
        feat_z = QgsFeature()
        feat_z.setGeometry(geom_p)
        feat_z.setAttributes([int(ad_forecast), "95% PI", len(fore_data)])
        pr_fore_zon.addFeatures([feat_z]); capa_fore_zon.updateExtents()
    else:
        print("  ⚠ Polígono de forecast vacío — revisar dirección del baseline")
else:
    print(f"  ⚠ Solo {len(fore_data)} transectas con PI — se necesitan ≥3")


# ═════════════════════════════════════════════════════════════════════════════
#  ESTILOS (9-bin DSAS)
# ═════════════════════════════════════════════════════════════════════════════

epr_vals = [f["EPR"] for f in capa_epr.getFeatures()]
lrr_vals = [f["LRR"] for f in capa_lrr.getFeatures()]
nsm_vals = [f["NSM"] for f in capa_nsm.getFeatures()]
sce_vals = [f["SCE"] for f in capa_sce.getFeatures()]
wlr_vals = [f["WLR"] for f in capa_wlr.getFeatures()]

capa_epr.setRenderer(renderer_dsas_9bins("EPR", epr_vals))
capa_lrr.setRenderer(renderer_dsas_9bins("LRR", lrr_vals))
capa_nsm.setRenderer(renderer_dsas_9bins("NSM", nsm_vals))
capa_wlr.setRenderer(renderer_dsas_9bins("WLR", wlr_vals))

# SCE siempre positivo → 5 bins azules
sce_max = max((v for v in sce_vals if v is not None), default=10)
sce_paso = sce_max / 5
colores_sce = ["#eff3ff","#bdd7e7","#6baed6","#2171b5","#08306b"]
rangos_sce  = [QgsRendererRange(i*sce_paso,(i+1)*sce_paso,
               QgsLineSymbol.createSimple({"color":colores_sce[i],"width":"0.9"}),
               f"{i*sce_paso:.1f}–{(i+1)*sce_paso:.1f} m") for i in range(5)]
capa_sce.setRenderer(QgsGraduatedSymbolRenderer("SCE", rangos_sce))

# Forecast
sym_lin = QgsLineSymbol.createSimple({"color":"#1a237e","width":"1.4","line_style":"solid"})
capa_fore_lin.renderer().setSymbol(sym_lin)
sym_zon = QgsFillSymbol.createSimple({"color":"70,130,180,70","style":"solid",
                                       "outline_color":"#1a237e","outline_width":"0.6"})
capa_fore_zon.renderer().setSymbol(sym_zon)

# ── Agregar al proyecto ──────────────────────────────────────────────────────
for capa in [capa_fore_zon, capa_fore_lin, capa_sce, capa_nsm,
             capa_wlr, capa_epr, capa_lrr]:
    QgsProject.instance().addMapLayer(capa)


# ═════════════════════════════════════════════════════════════════════════════
#  GUARDAR SHAPEFILES
# ═════════════════════════════════════════════════════════════════════════════

if CARPETA_SALIDA:
    mapeo = {"DSAS_EPR.shp":capa_epr, "DSAS_LRR.shp":capa_lrr,
             "DSAS_NSM.shp":capa_nsm, "DSAS_SCE.shp":capa_sce,
             "DSAS_WLR.shp":capa_wlr, "DSAS_Forecast_Lin.shp":capa_fore_lin,
             "DSAS_Forecast_Zona.shp":capa_fore_zon}
    print("\nGuardando shapefiles...")
    for nombre, capa in mapeo.items():
        ruta = os.path.join(CARPETA_SALIDA, nombre)
        err  = QgsVectorFileWriter.writeAsVectorFormat(
            capa, ruta, "UTF-8", baseline_capa.crs(), "ESRI Shapefile")
        estado = "✓" if err[0]==QgsVectorFileWriter.NoError else f"⚠ {err}"
        print(f"  {estado}  {nombre}")


# ═════════════════════════════════════════════════════════════════════════════
#  RESUMEN CONSOLA  (estadísticas DSAS — sin umbrales arbitrarios)
# ═════════════════════════════════════════════════════════════════════════════

lrr_validos  = [v for v in lrr_vals if v is not None]
n_lrr        = len(lrr_validos)
lrr_prom_res = sum(lrr_validos) / n_lrr if n_lrr else 0
lrr_sorted   = sorted(lrr_validos)
lrr_med      = (lrr_sorted[n_lrr//2-1]+lrr_sorted[n_lrr//2])/2 if n_lrr%2==0 else lrr_sorted[n_lrr//2] if n_lrr else 0
lrr_neg      = sum(1 for v in lrr_validos if v < 0)
lrr_pos      = sum(1 for v in lrr_validos if v >= 0)

epr_validos  = [v for v in epr_vals if v is not None]
epr_prom     = sum(epr_validos)/len(epr_validos) if epr_validos else 0

print(f"""
╔══════════════════════════════════════════════════════════════════╗
  RESULTADO — {total} transectas  |  {n_fechas} shorelines  |  {span:.1f} años
╠══════════════════════════════════════════════════════════════════╣
  EPR  promedio : {epr_prom:+.3f} m/año
  LRR  promedio : {lrr_prom_res:+.3f} m/año  |  mediana: {lrr_med:+.3f} m/año
       rango    : {min(lrr_validos):+.2f} a {max(lrr_validos):+.2f} m/año
╠══════════════════════════════════════════════════════════════════╣
  LRR < 0 (retroceso) : {lrr_neg:4d} transectas  ({lrr_neg/total*100:.1f}%)
  LRR ≥ 0 (avance)   : {lrr_pos:4d} transectas  ({lrr_pos/total*100:.1f}%)
╠══════════════════════════════════════════════════════════════════╣
  FORECAST {int(ad_forecast)} (+{AÑOS_FORECAST} a, PI 95%):
    Línea proyectada   → DSAS_Forecast_Lin
    Zona incertidumbre → DSAS_Forecast_Zona
╚══════════════════════════════════════════════════════════════════╝
""")
