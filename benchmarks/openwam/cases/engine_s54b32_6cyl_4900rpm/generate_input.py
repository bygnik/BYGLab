"""Generate the full 6-cylinder S54B32 OpenWAM input file."""
import math
import numpy as np

RPM = 4900
OUT = r"d:\Desktop\BYGLab\benchmarks\openwam\cases\engine_s54b32_6cyl_4900rpm\engine_s54b32_6cyl_4900rpm.txt"

lines = []
def L(text, comment=None):
    if comment:
        # guard: never allow a literal '>' inside the comment (breaks CleanLabelsX)
        assert '>' not in comment, f"stray > in comment: {comment}"
        lines.append(f"{text}\t<{comment}>")
    else:
        lines.append(text)

def blank():
    lines.append("")

# ---------- geometry / valve table precompute ----------
bore, stroke, CR = 0.087, 0.091, 11.5
area = math.pi * (bore/2)**2
Vs = area*stroke
Vc = Vs/(CR-1)
V_bdc = Vs+Vc
R, P0, T0 = 287.0, 1.05e5, 293.15
mass_seed = P0*V_bdc/(R*T0)

n_lev = 27
i = np.arange(n_lev)
lift = 0.012*0.5*(1-np.cos(2*np.pi*i/(n_lev-1)))
lift[0] = 1e-6
lift[-1] = 1e-6
lift_str = " ".join(f"{v:.6f}" for v in lift)

n_cd = 13
j = np.arange(n_cd)
lev = j*0.001
cd = np.minimum(0.64, 0.64*lev/0.004)
cd_str = " ".join(f"{v:.4f}" for v in cd)
torb_str = " ".join("0.0" for _ in range(n_cd))

# ============================================================
# Header
# ============================================================
L("2200", "OpenWAMVersion")
L("1", "Independent=true")
blank()

# ============================================================
# ReadGeneralData
# ============================================================
L("0.00001  2.0", "agincr(s) SimulationDuration=2.0 cycles (1.0 produced no INS.DAT: no cylinder completes a full local 720deg within 1 engine cycle given the firing offsets, and the output-flush never triggers)")
L("1.0     20.0", "AmbientPressure(bar) AmbientTemperature(degC)")
L("0       0", "tipocalculoespecies=0(nmCalculoSimple) tipogamma=0(nmGammaConstante)")
L("1", "hayBQ=1, EngineBlock=true")
L("0  0  0", "tipociclo=0(4T) tipomod=0(nmEstacionario) EGR=0")
L("1", "haycombustible=1")
L("1", "tipocombustible=1(gasoline/MEP), SpeciesNumber=4")
L("0.0  0.0  1.0", "CompAtmosfera: GasesQuemados=0 Gasolina=0 Aire=1.0")
blank()

# ============================================================
# ReadEngine -> LeeMotor
# ============================================================
L("0", "ACT flag = 0")
L("6", "NCilin = 6")
L(str(RPM), "FRegimen (RPM)")
L("1.05", "FPresionInicialRCA (bar, seed)")
L(f"{mass_seed:.6e}", "FMasaInicial (kg, seed mass, shared init for all 6 cylinders)")
L("0", "ImponerComposicionAE = 0")
L("0.0  0.0  1.0", "initial in-cylinder composition: GasesQuemados=0 Gasolina=0 Aire=1.0")
L("0", "TipoPresionAAE = 0")
L("1", "tipocombustion = 1 (gasoline/MEP)")
L("1.0", "FDosadoInicial, 1.0=stoichiometric")
L("0.98", "FRendimientoCombustion")
L("43000000", "FPoderCalorifico (J/kg)")
L("750", "FDensidadCombustible (kg/m3)")
L("1", "FNumTuboRendVol = 1 (cyl1 intake pipe, reporting reference only)")
L("450  550  420", "FTempInicial: Piston(K) Culata(K) Cylinder(K)")
L(f"{area:.7f}  {area:.7f}", "AreaPiston AreaCulata (m2)")
L("0.01  150  2700  900", "FParedPiston (inert, CalculoTempPared=2)")
L("0.01  150  2700  900", "FParedCulata (inert)")
L("0.01  150  2700  900", "FParedCilindro (inert)")
L("1.0  1.0  100000  90", "FAjusteTranCalAdm FAjusteTranCalEsc FParPotMax(inert) FTempRefrigerante(inert)")
L("2", "CalculoTempPared = 2 (nmTempFija, fixed wall temp)")
L("2.28  0.00324  0.0", "Woschni cw1 cw2 xpe(dead field)")
L("0.139", "FGeom.Biela (m)")
L("0.091", "FGeom.Carrera (m)")
L("0.087", "FGeom.Diametro (m)")
L("11.5", "FGeom.RelaCompresion (real S54B32 spec)")
L("0.0", "FGeom.DiametroBowl = 0")
L("0.0", "FGeom.AlturaBowl")
L("0.02", "FGeom.DistanciaValvulas")
L("0.0", "FGeom.AreaBlowBy = 0")
L("0.7", "FGeom.CDBlowBy (inert)")
L("0.0", "FGeom.Excentricidad (mm) = 0")
L("0.021", "FGeom.DiametroBulon (inert)")
L("0.03", "FGeom.AlturaCoronaPiston (inert)")
L("0.45", "FGeom.MasaBiela (inert)")
L("0.35", "FGeom.MasaPistonSegmentosBulon (inert)")
L("2.1e11", "FGeom.ModuloElasticidad (inert)")
L("0.0", "FGeom.CoefDeformaciones = 0")
L("0  0  0  0", "FPerMec.Coef0..3")
L("1", "FNumeroLeyesQuemado = 1")
L("0  0  0", "ma mf n (irrelevant, size==1)")
L("1", "FNWiebes = 1")
L("2.5  6.9  1.0  45  15", "Wiebe: m C Beta IncAlpha(deg) Alpha0(deg BTDC) - shared by all 6 cylinders")
L("0", "FTipoDatosIny = 0")
# --- multi-cylinder firing order block (NCilin=6, tipodesfa=1=nmImpuesto) ---
L("1", "tipodesfa = 1 (nmImpuesto, firing-order permutation)")
L("1", "firing position 0: cylinder 1 (FDesfase[0]=0)")
L("5", "firing position 1: cylinder 5 (FDesfase[4]=120)")
L("3", "firing position 2: cylinder 3 (FDesfase[2]=240)")
L("6", "firing position 3: cylinder 6 (FDesfase[5]=360)")
L("2", "firing position 4: cylinder 2 (FDesfase[1]=480)")
L("4", "firing position 5: cylinder 4 (FDesfase[3]=600)")
L("0", "global engine-level controllers count = 0")
for c in range(1, 7):
    L("0", f"cyl {c} per-cylinder controllers count = 0")
blank()

# ============================================================
# ReadPipes: 15 pipes
#   1-6:   cyl1-6 intake runners (ITB, atmosphere -> intake valve)
#   7,8,9: cyl1,3,2 exhaust runners -> branch A (group firing 240deg apart: 1,3,2)
#   10:    secondary A -> branch C
#   11,12,13: cyl5,6,4 exhaust runners -> branch B (group: 5,6,4)
#   14:    secondary B -> branch C
#   15:    tailpipe branch C -> atmosphere
# ============================================================
NumberOfPipes = 15
L(str(NumberOfPipes), "NumberOfPipes")

def pipe_general(nodoizq, nododer, comment):
    L(f"1  {nodoizq}  {nododer}  1  1" if False else f"{nodoizq}  {nododer}  1  1", comment)
    L("0.0", "FFriccion")
    L("20.0  20.0  1.0  0.0", "FTIniParedTub FTini(degC) FPini=1.0bar FVelMedia=0")
    L("1  0.0  0.0", "TipTC=1 FCoefAjusTC=0 FCoefAjusFric=0")
    L("0.0  0.0  1.0", "FComposicionInicial: air")
    L("0.01  2", "FMallado=0.01m FTctpt=2")
    L("2", "metodo[0]=2(nmTVD)")
    L("0.4", "FCourant")

def pipe_geom(d0, length, d1, comment):
    L(f"{d0}", comment)
    L(f"{length}  {d1}", "FLTramo[1] FDExtTramo[1]")

# intake runners: cyl 1..6, node pairs (1,2)(3,4)(5,6)(7,8)(9,10)(11,12)
intake_nodes = [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12)]
for idx, (nl, nr) in enumerate(intake_nodes, start=1):
    pipe_general(nl, nr, f"Pipe{idx}: cyl{idx} intake runner (ITB), atm(node{nl}) to intake valve(node{nr})")
    pipe_geom(0.045, 0.3, 0.045, "FDExtTramo[0]=0.045m")

# exhaust primary runners group A: cyl1(node13->14) cyl3(node15->14) cyl2(node16->14)
pipe_general(13, 14, "Pipe7: cyl1 exhaust runner to branch A")
pipe_geom(0.04, 0.35, 0.04, "FDExtTramo[0]=0.04m")
pipe_general(15, 14, "Pipe8: cyl3 exhaust runner to branch A")
pipe_geom(0.04, 0.35, 0.04, "FDExtTramo[0]=0.04m")
pipe_general(16, 14, "Pipe9: cyl2 exhaust runner to branch A")
pipe_geom(0.04, 0.35, 0.04, "FDExtTramo[0]=0.04m")

# secondary A: branch A(14) -> branch C(17), larger diameter
pipe_general(14, 17, "Pipe10: secondary A, branch A(node14) to branch C(node17)")
pipe_geom(0.055, 0.3, 0.055, "FDExtTramo[0]=0.055m")

# exhaust primary runners group B: cyl5(node18->19) cyl6(node20->19) cyl4(node21->19)
pipe_general(18, 19, "Pipe11: cyl5 exhaust runner to branch B")
pipe_geom(0.04, 0.35, 0.04, "FDExtTramo[0]=0.04m")
pipe_general(20, 19, "Pipe12: cyl6 exhaust runner to branch B")
pipe_geom(0.04, 0.35, 0.04, "FDExtTramo[0]=0.04m")
pipe_general(21, 19, "Pipe13: cyl4 exhaust runner to branch B")
pipe_geom(0.04, 0.35, 0.04, "FDExtTramo[0]=0.04m")

# secondary B: branch B(19) -> branch C(17)
pipe_general(19, 17, "Pipe14: secondary B, branch B(node19) to branch C(node17)")
pipe_geom(0.055, 0.3, 0.055, "FDExtTramo[0]=0.055m")

# tailpipe: branch C(17) -> atmosphere(22)
pipe_general(17, 22, "Pipe15: tailpipe, branch C(node17) to atmosphere(node22)")
pipe_geom(0.065, 0.5, 0.065, "FDExtTramo[0]=0.065m")

blank()

# ============================================================
# ReadValves: 12 valves (6 intake identical, 6 exhaust identical)
# valve index 1-6 = intake cyl1-6, valve index 7-12 = exhaust cyl1-6
# ============================================================
L("12", "NumberOfValves")

def valve(diam, angulo_apertura, comment):
    L("1", f"{comment}: type=1 (TValvula4T)")
    L(f"{diam}", "FDiametro (m)")
    L("27", "NumLev")
    L("10.0", "FIncrAng (deg)")
    L(f"{angulo_apertura}", "FAnguloApertura (deg)")
    L("0", "FDiametroRef = 0, no Cd rescale")
    L("0", "FCoefTorbMedio (dead field)")
    L(lift_str, "FLevantamiento, 27 values")
    L("13", "NumCD")
    L("0.001", "FIncrLev (m)")
    L(cd_str, "FDatosCDEntrada, 13 values")
    L(cd_str, "FDatosCDSalida, 13 values")
    L(torb_str, "FDatosTorbellino, 13 values")
    L("1", "ControlRegimen = 1, cylinder-bound")
    L("1.0", "FRelacionVelocidades (inert)")
    L("0", "VVT controllers count = 0")

for c in range(1, 7):
    valve(0.0495, 340.0, f"Valve{c}: intake cyl{c} (IVO 340deg)")
for c in range(1, 7):
    valve(0.0431, 140.0, f"Valve{6+c}: exhaust cyl{c} (EVO 140deg)")
blank()

# ============================================================
# ReadPlenums / ReadCompressors
# ============================================================
L("0", "NumberOfPlenums")
L("0  0  0", "WAMer-only counters")
blank()
L("0", "NumberOfCompressors")
blank()

# ============================================================
# ReadConnections: 22 connections
#   1:  atm (Pipe1 left, cyl1 intake)
#   2:  intake valve cyl1 (Pipe1 right)               quevalv=1
#   3:  atm (Pipe2 left, cyl2 intake)
#   4:  intake valve cyl2                              quevalv=2
#   5:  atm (Pipe3 left, cyl3 intake)
#   6:  intake valve cyl3                              quevalv=3
#   7:  atm (Pipe4 left, cyl4 intake)
#   8:  intake valve cyl4                              quevalv=4
#   9:  atm (Pipe5 left, cyl5 intake)
#   10: intake valve cyl5                               quevalv=5
#   11: atm (Pipe6 left, cyl6 intake)
#   12: intake valve cyl6                               quevalv=6
#   13: exhaust valve cyl1 (Pipe7 left)                 quevalv=7
#   14: branch A (Pipe7,8,9 right ends + Pipe10 left)   TipoCC=12
#   15: exhaust valve cyl3 (Pipe8 left)                 quevalv=9
#   16: exhaust valve cyl2 (Pipe9 left)                 quevalv=8
#   17: branch C (Pipe10,14 right ends + Pipe15 left)   TipoCC=12
#   18: exhaust valve cyl5 (Pipe11 left)                quevalv=11
#   19: branch B (Pipe11,12,13 right ends + Pipe14 left) TipoCC=12
#   20: exhaust valve cyl6 (Pipe12 left)                quevalv=12
#   21: exhaust valve cyl4 (Pipe13 left)                quevalv=10
#   22: atm (Pipe15 right, tailpipe)
# ============================================================
NumberOfConnections = 22
L(str(NumberOfConnections), "NumberOfConnections")
L("0 0 0 0 0 0 0 0 0", "WAMer-only counters")

def atm(node):
    L("0", f"BC#{node}: TipoCC=0 (nmOpenEndAtmosphere)")
    L("1.0", "FPerdidaExtremo")

def intake_valve(node, cyl, quevalv):
    L("7", f"BC#{node}: TipoCC=7 (nmIntakeValve) cyl{cyl}")
    L(f"1  {cyl}", "numid(dead) FNumeroCilindro")
    L(f"{quevalv}", f"quevalv={quevalv} (intake valve def for cyl{cyl})")

def exhaust_valve(node, cyl, quevalv):
    L("8", f"BC#{node}: TipoCC=8 (nmExhaustValve) cyl{cyl}")
    L(f"1  {cyl}", "numid(dead) FNumeroCilindro")
    L(f"{quevalv}", f"quevalv={quevalv} (exhaust valve def for cyl{cyl})")

def branch(node, label):
    L("12", f"BC#{node}: TipoCC=12 (nmBranch) {label}")

atm(1)
intake_valve(2, 1, 1)
atm(3)
intake_valve(4, 2, 2)
atm(5)
intake_valve(6, 3, 3)
atm(7)
intake_valve(8, 4, 4)
atm(9)
intake_valve(10, 5, 5)
atm(11)
intake_valve(12, 6, 6)
exhaust_valve(13, 1, 7)
branch(14, "branch A: cyl1+cyl3+cyl2 exhaust runners + secondary A")
exhaust_valve(15, 3, 9)
exhaust_valve(16, 2, 8)
branch(17, "branch C: secondary A + secondary B + tailpipe")
exhaust_valve(18, 5, 11)
branch(19, "branch B: cyl5+cyl6+cyl4 exhaust runners + secondary B")
exhaust_valve(20, 6, 12)
exhaust_valve(21, 4, 10)
atm(22)
blank()

# ============================================================
# Axis / Sensors / Controllers
# ============================================================
L("0", "NumberOfAxis")
blank()
L("0", "NumberOfSensors")
blank()
L("0", "NumberOfControllers")
blank()

# ============================================================
# Output: Average results (must set NumEnginesAvg>0 to avoid AvgEngine NULL deref bug)
# ============================================================
L("2", "FTypeOfInsResults=2 (nmAllCyclesConcatenated) - this field actually governs INSTANTANEOUS results retention despite living in the AVG-results read function; 0=nmLastCycle silently keeps only the final cycle's data, truncating the file")
L("0", "NumPipesAvg")
L("0", "WAMer dummy int after NumPipesAvg")
L("0", "NumCylindersAvg")
L("1", "NumEnginesAvg = 1 (required workaround)")
L("2", "engine avg nvars")
L("1 12", "engine avg var IDs: 1=ParEfectivo(torque) 12=Potencia(power)")
L("0", "NumPlenumsAvg")
L("0", "NumAxisAvg")
L("0", "NumCompressorAvg")
L("0", "NumTurbineAvg")
L("0", "NumValvesAvg")
L("0", "NumRootsAvg")
L("0", "NumVenturisAvg")
L("0", "NumConnectionsAvg")
L("0", "NumSensorAvg")
L("0", "NumControllersAvg")
blank()

# ============================================================
# Output: Instantaneous results -- all 6 cylinders' P/T/V traces
# ============================================================
L("6", "NumCylindersIns = 6 (all cylinders)")
for c in range(1, 7):
    L(f"{c}", f"CylinderID = {c}")
    L("3", "FNumVarIns = 3")
    L("0 1 11", "var IDs: 0=Pressure 1=Temperature 11=Volumen")
L("0", "NumPlenumsIns")
L("0", "NumPipesIns")
L("0", "WAMer dummy int after NumPipesIns")
L("0", "NumVenturisIns")
L("0", "NumValvesIns")
L("0", "NumTurboIns")
L("0", "NumCompressorIns")
L("0", "NumTurbineIns")
L("0", "NumRootsIns")
L("0", "NumConnectionsIns")
L("0", "NumWasteGateIns")
L("0", "NumReedIns")
L("0", "NumSensorIns")
L("0", "NumControllersIns")
blank()

# ============================================================
# Space-time results: skip (0), rely on per-cylinder INS.DAT channel
# ============================================================
L("0", "FNumMagnitudesEspTemp: 0, skip per-cell space-time dump for this case")
blank()

L("0", "dll=0, ReadDataDLL skipped")

with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")

print(f"Wrote {len(lines)} lines to {OUT}")
