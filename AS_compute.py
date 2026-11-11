#%% Preamble
from IPython.display import display, Math, Latex
from sympy import * #load symbolic library
from sympy.polys.matrices import DomainMatrix
import os
import pickle

# =========================
# FOLDERS/PKLs
# =========================
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

# =========================
# SYMBOLS
# =========================
q  = Symbol('q')
s  = Symbol('s')
x=symbols('x0:4')
y=symbols('y0:4')

def totensor(x):
    return (x//4,x%4)

def condense(R):
    C = zeros(4)
    for j in range(16):
        for i in range(16):
            xind,yind = totensor(i)
            entry = factor(R[i,j])
            C[totensor(j)] += expand(numer(entry))/expand(denom(entry))*x[xind]*y[yind]
    return C

# =========================
# WORKER FUNCTIONS
# =========================
def sparseM(x):
    return  ImmutableSparseMatrix(x)

def KP(*factors):
    prod = Matrix([[1]])
    for x in factors:
        prod = kronecker_product(prod,x)
    return  sparseM(prod)

def get_LT(M, x):
    def leading_term(expr):
        terms = Add.make_args(expand(expr))
        degrees = [t.as_powers_dict().get(x, 0) for t in terms]
        max_degree = max(degrees)

        return Add(*[
            t for t, d in zip(terms, degrees)
            if d == max_degree
        ])

    return M.applyfunc(leading_term)

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
    assert M == m*eye(*shape(M))
    return m

double = RT['l']*RT['r']
doublei = (double.to_DM().inv()).to_Matrix()
cupcap = RT['cl']*RT['er']
crL = filler(0,1,'cr')
elL = filler(0,3,'el')
erU = filler(0,3,'er')
crU = filler(1,1,'cr')
crR = crU
elR = filler(1,3,'el')

HR = matrix_scalar(elR*KP(double, eye(4))*crR)
# HHR = matrix_scalar(elR*KP(double2, eye(4))*crR)

rels = Matrix([[0,1,HR], [1,0,1], [q**2*HR,1,0]])
relsi = rels.inv()

name_to_index = {'id':0, 'cupcap':1, 'double':2}
name_to_diagram = {'id':eye(16), 'cupcap':cupcap, 'double':double}
basis = ['id', 'cupcap', 'double']

def tangle_2_cmpnts(M):
    ML = elR*KP(M, eye(4))*crR
    ml = matrix_scalar(ML)

    MU = erU * KP(M,eye(4)) * crU
    mu = matrix_scalar(MU)

    MD = elR*KP(doublei*M, eye(4))*crR
    md = matrix_scalar(MD)

    v = Matrix([ml, mu, md])
    T = simplify(relsi*v)
    check = T[0]*eye(16)+ T[1]*cupcap + T[2]*double
    assert simplify(expand(M-check)) == zeros(*shape(M))
    return T


fname = os.path.join(cwd, 'pkls',"trefoil_6twists.pkl")
# tangle = "0cl, 1cr, 3li, 2u, 4d, 3ri, 3li, 2u, 4d, 3ri, 3li, 2u, 4d, 3ri, 1el, 0er, 0r, 0l, 0r, 0l, 0r, 0l "
# bottom_width = 2
trefoiltwist = disk_load(fname).to_DM()

fname = os.path.join(cwd, 'pkls',"pclasp.pkl")
# tangle = "1cl, 0u, 2d, 1er"
# bottom_width = 2
pclasp = disk_load(fname).to_DM()

fname = os.path.join(cwd, 'pkls',"mclasp.pkl")
# tangle = "1cl, 0ui, 2di, 1er"
# bottom_width = 2
mclasp = disk_load(fname).to_DM()

fname = os.path.join(cwd, 'pkls',"AS_top.pkl")
# tangle = '1r,0er'
# bottom_width = 3
AS_top = disk_load(fname)

fname = os.path.join(cwd, 'pkls',"AS_bottom.pkl")
# tangle = '1cl, 2cr, 1r, 0r,0l,1l,2el'
# bottom_width = 1
AS_bottom = disk_load(fname)

bottomcr = filler(1,2,'cr')
topel = filler(1,4,'el')

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


final_tangle = zeros(3,1)
for a in basis:
    M = AS_top*KP(eye(4),name_to_diagram[a])*AS_bottom
    assert simplify(M[0,0]*eye(4)-M) == zeros(4)
    final_tangle[name_to_index[a]] = simplify(M[0,0])

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

'''
Print results:
'''
cc = condense(cupcap)
display(Math(rf"cupcap=({latex(factor(cc[0,0]))}){latex(simplify(cc/cc[0,0]))}"))

doub = condense(double)
display(Math(rf"double={latex(doub)}"))

collectttp = ttp_components.applyfunc(lambda x: collect(x,s))
display(Math(rf"TT^+={latex(ttp_components)}"))

collectttm = ttm_components.applyfunc(lambda x: collect(x,s))
display(Math(rf"TT^-={latex(ttm_components)}"))

for a in basis:
    for b in basis:
        X=latex(join_AB(eye(3)[:,name_to_index[a]],eye(3)[:,name_to_index[b]]).applyfunc(factor).applyfunc(lambda x: collect(x,s)))
        display(Math(rf"\text{{{a}}}\boxtimes\text{{{b}}}={X}"))

for a in basis:
    X=latex(collect(final_tangle[name_to_index[a]],s))
    display(Math(rf"\mathrm{{AS}}^*(\text{{{a}}})={X}"))

X = latex(collect(AS_link1,s))
display(Math(rf"\mathrm{{LG}}_{{\mathrm{{AS}}(1)}}(s,q)={X}"))

X = latex(factor(AS_link1.subs(q,1)))
display(Math(rf"\mathrm{{LG}}_{{\mathrm{{AS}}(1)}}(s,1)={X}"))

T = get_LT(get_LT(fatgraph,s),q)
X = latex(join_AB(Matrix([x[0],x[1],x[2]]),T).applyfunc(lambda expr: collect(expr, [x[0],x[1],x[2]])))
Y = latex(Matrix([x[0],x[1],x[2]]))
display(Math(rf"{Y}\boxtimes T\sim {X}"))
