#%% Preamble
from sympy import * #load symbolic library
from sympy.polys.matrices import DomainMatrix
import os
import pickle

cwd = os.getcwd()
pkl_dir = os.path.join(cwd, "pkls")

def disk_load(fname):
    if os.path.exists(fname):
        with open(fname, "rb") as f:
            return pickle.load(f)
    return None

def disk_save(fname, obj):
    with open(fname, "wb") as f:
        pickle.dump(obj, f)

q  = Symbol('q')
s  = Symbol('s')

def sparseM(x):
    return  ImmutableSparseMatrix(x)

def KP(*factors):
    prod = Matrix([[1]])
    for x in factors:
        prod = kronecker_product(prod,x)
    return  sparseM(prod)

def saved_computation(name, function, args=None):
    if args is None:
        args = ()
    elif not isinstance(args, tuple):
        args = (args,)
    args_str = ",".join(str(a) for a in args)
    fname = os.path.join(pkl_dir, f"{name}({args_str}).pkl")

    X = disk_load(fname)
    if X is not None:
        return X
    X = function(*args)
    disk_save(fname, X)
    return X

dim = 4
fname = os.path.join(cwd,'pkls', 'LG_RT_data.pkl')
RT = disk_load(fname)

def filler(left_strands, bottom_width, elementary_tangle):
    fname = os.path.join(cwd,pkl_dir,f"RT_{left_strands}_{bottom_width}_{elementary_tangle}.pkl")
    slice_diagram = disk_load(fname)
    if slice_diagram is not None:
        # return slice_diagram.to_DM()
        return slice_diagram
    right_strands = bottom_width - left_strands - 2
    if elementary_tangle in {"cl", "cr"}:
        right_strands = bottom_width - left_strands
    if elementary_tangle in {'id'}:
        right_strands = bottom_width - left_strands - 1
    slice_diagram = KP(eye(dim**left_strands), RT[elementary_tangle], eye(dim**right_strands))
    disk_save(fname, slice_diagram)
    # slice_diagram = slice_diagram.to_DM()
    return slice_diagram

def matrix_scalar(M):
    M = simplify(M)
    m = M[0,0]
    if M == m*eye(*shape(M)):
        return m
    print("Matrix is not scalar")

double = RT['l']*RT['r']
double2 = (double.to_DM()**2).to_Matrix()
cupcap = RT['cl']*RT['er']
crL = filler(0,1,'cr')
elL = filler(0,3,'el')
erU = filler(0,3,'er')
crU = filler(1,1,'cr')
crR = crU
elR = filler(1,3,'el')

HR = matrix_scalar(elR*KP(double, eye(4))*crR)
HHR = matrix_scalar(elR*KP(double2, eye(4))*crR)

rels = Matrix([[0,1,HR], [1,0,1], [HR,1,HHR]])
relsi = rels.inv()

name_to_index = {'id':0, 'cupcap':1, 'double':2}
name_to_diagram = {'id':eye(16), 'cupcap':cupcap, 'double':double}
basis = ['id', 'cupcap', 'double']

def tangle_2_cmpnts(M):
    ML = elL*KP(eye(4),M)*crL
    ml = matrix_scalar(ML)

    MU = erU * KP(M,eye(4)) * crU
    mu = matrix_scalar(MU)

    MD = elL*KP(eye(4),double*M)*crL
    md = matrix_scalar(MD)

    v = Matrix([ml, mu, md])
    T = simplify(relsi*v)
    check = T[0]*eye(16)+ T[1]*cupcap + T[2]*double
    if simplify(expand(M-check)) == zeros(*shape(M)):
        return T
    print("Error")
    return


fname = os.path.join(cwd, 'pkls',"trefoil_6twists.pkl")
trefoiltwist = disk_load(fname).to_DM()

fname = os.path.join(cwd, 'pkls',"pclasp.pkl")
pclasp = disk_load(fname).to_DM()

fname = os.path.join(cwd, 'pkls',"mclasp.pkl")
mclasp = disk_load(fname).to_DM()

fname = os.path.join(cwd, 'pkls',"cr_at_left_1_bottom_2.pkl")
bottomcr = disk_load(fname)
fname = os.path.join(cwd, 'pkls',"el_at_left_1_bottom_4.pkl")
topel = disk_load(fname)

fname = os.path.join(cwd, 'pkls',"AS_top.pkl")
AS_top = disk_load(fname)

fname = os.path.join(cwd, 'pkls',"AS_bottom.pkl")
AS_bottom = disk_load(fname)

# ttp = (pclasp*trefoiltwist).to_Matrix()
# ttm = (mclasp*trefoiltwist).to_Matrix()

# ttp_components = tangle_2_cmpnts(ttp)
fname = os.path.join(cwd, 'pkls',"ttp_components.pkl")
# disk_save(fname, ttp_components)
ttp_components = disk_load(fname)

# ttm_components = tangle_2_cmpnts(ttm)
fname = os.path.join(cwd,pkl_dir,"ttm_components.pkl")
# disk_save(fname, ttm_components)
ttm_components = disk_load(fname)

def compute_fat_graph_comps():
    fat_graph_comps = {}
    for a in basis:
        for b in basis:
            M = topel*KP(name_to_diagram[a],name_to_diagram[b])*bottomcr
            fat_graph_comps[(a,b)] = expand(tangle_2_cmpnts(M))
    return fat_graph_comps

fat_graph_comps = saved_computation('fat_graph_comps', compute_fat_graph_comps)

def join_AB(A,B):
    components = zeros(3,1)
    for a in basis:
        for b in basis:
            components += fat_graph_comps[a,b]*A[name_to_index[a]]*B[name_to_index[b]]
    return expand(components)

fatgraph = join_AB(ttp_components,ttm_components)


final_tangle_id = (s - 1)*(s + 1)*(q*s - 1)*(q*s + 1)*(q**4*s**8 + q**4*s**4 + q**2*s**4 + 1)/(q**6*s**6)
final_tangle_cupcap = (s - 1)**2*(s + 1)**2*(q*s - 1)**2*(q*s + 1)**2/(q**4*s**4)
final_tangle_double = (s - 1)*(s + 1)*(q*s - 1)*(q*s + 1)*(3*q**6*s**8 - 4*q**6*s**6 + 3*q**6*s**4 + q**4*s**8 - 4*q**4*s**6 + 7*q**4*s**4 - 4*q**4*s**2 + 2*q**2*s**4 - 4*q**2*s**2 + 3*q**2 + 1)/(q**8*s**6)

name_to_final_tangle = {'id':final_tangle_id, 'cupcap':final_tangle_cupcap, 'double':final_tangle_double}
final_tangle = zeros(3,1)
for a in basis:
    final_tangle[name_to_index[a]] = simplify(name_to_final_tangle[a])

def AS_link(n):
    components = fatgraph
    for i in range(n-1):
        components = join_AB(components,fatgraph)
    output = sum(components[name_to_index[a]]*final_tangle[name_to_index[a]] for a in basis)
    return expand(output)

AS_link1 = saved_computation('AS_link', AS_link, 1)
AS_link2 = saved_computation('AS_link', AS_link, 2)
AS_link3 = saved_computation('AS_link', AS_link, 3)
AS_link4 = saved_computation('AS_link', AS_link, 4)


