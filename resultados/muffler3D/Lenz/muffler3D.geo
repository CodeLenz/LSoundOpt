// =============================================================================
// SCRIPT GMSH - MUFFLER CAIXA 3D (100% HEXAEDROS)
// ============================================================================

lc = 0.01;

// Dimensões Y (Altura)
H_tube   = 0.10; 
H_cavity = 0.20;

// Dimensões Z (Profundidade) - Seção Quadrada/Cúbica
Z_tube   = 0.10; // Tubo quadrado (0.10 x 0.10)
Z_cavity = 0.20;

// Câmara em caixa envolvendo o tubo
// Dimensões X (Comprimento)
L_in      = 0.4;  
L_chamber = 0.6;
L_gap     = 0.2;  
L_end     = 0.2;  

// Coordenadas X
X0 = 0.0;
X1 = L_in;
X2 = L_in + L_chamber;
X3 = X2 + L_gap;     
X4 = X3 + L_end;

// Coordenadas Z
Z0 = 0.0;
Z1 = Z_cavity;
Z2 = Z_cavity + Z_tube;
Z3 = Z_cavity + Z_tube + Z_cavity;

// Discretização (Nós por aresta)
n_in      = Ceil((X1 - X0) / lc) + 1;
n_chamber = Ceil((X2 - X1) / lc) + 1;
n_gap     = Ceil((X3 - X2) / lc) + 1;
n_end     = Ceil((X4 - X3) / lc) + 1;

n_z_cav   = Ceil(Z_cavity / lc) + 1;
n_z_tube  = Ceil(Z_tube / lc) + 1;

n_y_cav   = Ceil(H_cavity / lc) + 1;
n_y_tube  = Ceil(H_tube / lc) + 1;

// =============================================================================
// 1. BASE 2D NO PLANO X-Z (Criada na cota Y = H_cavity)
// =============================================================================
Y_base = H_cavity;

// PONTOS (Grid no plano horizontal)
p1 = newp; Point(p1) = {X0, Y_base, Z1, lc};
p4 = newp; Point(p4) = {X0, Y_base, Z2, lc};

p5 = newp; Point(p5) = {X1, Y_base, Z0, lc};
p2 = newp; Point(p2) = {X1, Y_base, Z1, lc};
p3 = newp; Point(p3) = {X1, Y_base, Z2, lc};
p9 = newp; Point(p9) = {X1, Y_base, Z3, lc};

p6  = newp; Point(p6)  = {X2, Y_base, Z0, lc};
p7  = newp; Point(p7)  = {X2, Y_base, Z1, lc};
p8  = newp; Point(p8)  = {X2, Y_base, Z2, lc};
p10 = newp; Point(p10) = {X2, Y_base, Z3, lc};

p11 = newp; Point(p11) = {X3, Y_base, Z1, lc};
p12 = newp; Point(p12) = {X3, Y_base, Z2, lc};

p13 = newp; Point(p13) = {X4, Y_base, Z1, lc};
p14 = newp; Point(p14) = {X4, Y_base, Z2, lc};

// LINHAS (X-dir)
lx_in_1 = newl; Line(lx_in_1) = {p1, p2};
lx_in_2 = newl; Line(lx_in_2) = {p4, p3};
lx_cb   = newl; Line(lx_cb)   = {p5, p6};
lx_cc1  = newl; Line(lx_cc1)  = {p2, p7};
lx_cc2  = newl; Line(lx_cc2)  = {p3, p8};
lx_cf   = newl; Line(lx_cf)   = {p9, p10};
lx_o1_1 = newl; Line(lx_o1_1) = {p7, p11};
lx_o1_2 = newl; Line(lx_o1_2) = {p8, p12};
lx_o2_1 = newl; Line(lx_o2_1) = {p11, p13};
lx_o2_2 = newl; Line(lx_o2_2) = {p12, p14};

// LINHAS (Z-dir)
lz_in = newl; Line(lz_in) = {p1, p4};
lz_c1 = newl; Line(lz_c1) = {p5, p2};
lz_c2 = newl; Line(lz_c2) = {p2, p3};
lz_c3 = newl; Line(lz_c3) = {p3, p9};
lz_c4 = newl; Line(lz_c4) = {p6, p7};
lz_c5 = newl; Line(lz_c5) = {p7, p8};
lz_c6 = newl; Line(lz_c6) = {p8, p10};
lz_o1 = newl; Line(lz_o1) = {p11, p12};
lz_o2 = newl; Line(lz_o2) = {p13, p14};

// LOOPS E SUPERFÍCIES (Planta Baixa)
cl_in = newll; Curve Loop(cl_in) = {lx_in_1, lz_c2, -lx_in_2, -lz_in};
s_in = news; Plane Surface(s_in) = {cl_in};

cl_cb = newll; Curve Loop(cl_cb) = {lx_cb, lz_c4, -lx_cc1, -lz_c1};
s_cb = news; Plane Surface(s_cb) = {cl_cb};

cl_cc = newll; Curve Loop(cl_cc) = {lx_cc1, lz_c5, -lx_cc2, -lz_c2};
s_cc = news; Plane Surface(s_cc) = {cl_cc};

cl_cf = newll; Curve Loop(cl_cf) = {lx_cc2, lz_c6, -lx_cf, -lz_c3};
s_cf = news; Plane Surface(s_cf) = {cl_cf};

cl_o1 = newll; Curve Loop(cl_o1) = {lx_o1_1, lz_o1, -lx_o1_2, -lz_c5};
s_o1 = news; Plane Surface(s_o1) = {cl_o1};

cl_o2 = newll; Curve Loop(cl_o2) = {lx_o2_1, lz_o2, -lx_o2_2, -lz_o1};
s_o2 = news; Plane Surface(s_o2) = {cl_o2};

// MALHA TRANSFINITA 2D (Estruturando o piso)
Transfinite Curve {lx_in_1, lx_in_2} = n_in;
Transfinite Curve {lx_cb, lx_cc1, lx_cc2, lx_cf} = n_chamber;
Transfinite Curve {lx_o1_1, lx_o1_2} = n_gap;
Transfinite Curve {lx_o2_1, lx_o2_2} = n_end;

Transfinite Curve {lz_c1, lz_c4, lz_c3, lz_c6} = n_z_cav;
Transfinite Curve {lz_in, lz_c2, lz_c5, lz_o1, lz_o2} = n_z_tube;

Transfinite Surface {s_in} = {p1, p2, p3, p4};
Transfinite Surface {s_cb} = {p5, p6, p7, p2};
Transfinite Surface {s_cc} = {p2, p7, p8, p3};
Transfinite Surface {s_cf} = {p3, p8, p10, p9};
Transfinite Surface {s_o1} = {p7, p11, p12, p8};
Transfinite Surface {s_o2} = {p11, p13, p14, p12};

Recombine Surface {s_in, s_cb, s_cc, s_cf, s_o1, s_o2};

// =============================================================================
// 2. EXTRUSÃO 3D EM Y (Stacking dos volumes)
// =============================================================================

// Extrusão para BAIXO (Layer 0: Fundo da Câmara)
e_bot_cb = Extrude {0, -H_cavity, 0} { Surface{s_cb}; Layers{n_y_cav}; Recombine; };
e_bot_cc = Extrude {0, -H_cavity, 0} { Surface{s_cc}; Layers{n_y_cav}; Recombine; };
e_bot_cf = Extrude {0, -H_cavity, 0} { Surface{s_cf}; Layers{n_y_cav}; Recombine; };

// Extrusão para CIMA (Layer 1: Caminho principal do Tubo)
e_mid_in = Extrude {0, H_tube, 0} { Surface{s_in}; Layers{n_y_tube}; Recombine; };
e_mid_cb = Extrude {0, H_tube, 0} { Surface{s_cb}; Layers{n_y_tube}; Recombine; };
e_mid_cc = Extrude {0, H_tube, 0} { Surface{s_cc}; Layers{n_y_tube}; Recombine; };
e_mid_cf = Extrude {0, H_tube, 0} { Surface{s_cf}; Layers{n_y_tube}; Recombine; };
e_mid_o1 = Extrude {0, H_tube, 0} { Surface{s_o1}; Layers{n_y_tube}; Recombine; };
e_mid_o2 = Extrude {0, H_tube, 0} { Surface{s_o2}; Layers{n_y_tube}; Recombine; };

// Extrusão para CIMA (Layer 2: Topo da Câmara)
// Usamos as faces do teto (índice 0) geradas na layer anterior
e_top_cb = Extrude {0, H_cavity, 0} { Surface{e_mid_cb[0]}; Layers{n_y_cav}; Recombine; };
e_top_cc = Extrude {0, H_cavity, 0} { Surface{e_mid_cc[0]}; Layers{n_y_cav}; Recombine; };
e_top_cf = Extrude {0, H_cavity, 0} { Surface{e_mid_cf[0]}; Layers{n_y_cav}; Recombine; };

// =============================================================================
// 3. COERÊNCIA E GRUPOS FÍSICOS (Imunes a deleção de IDs)
// =============================================================================

// Forçar a fusão topológica de faces coincidentes ANTES de criar os grupos
Coherence; 

// Volumes separados corretamente para evitar sobreposição na matriz de rigidez
Physical Volume("Material,Ar,1,1.21,343.0,0.0") = {
    e_bot_cb[1], e_bot_cc[1], e_bot_cf[1],
    e_mid_cb[1], e_mid_cf[1], 
    e_top_cb[1], e_top_cc[1], e_top_cf[1],
    e_mid_in[1], e_mid_cc[1], e_mid_o1[1], e_mid_o2[1]
};

Physical Volume("Fixed,0.0") = {
    e_mid_in[1], e_mid_cc[1], e_mid_o1[1], e_mid_o2[1]
};

// Tolerância para o BoundingBox
eps = 1e-4;

// Inlet (X=X0)
Physical Surface("Vn,1E-3,0.0,0.0") = Surface In BoundingBox{
    X0-eps, Y_base-eps, Z1-eps, 
    X0+eps, Y_base+H_tube+eps, Z2+eps
};

// Target / Outlet 1 / Gap Interface (X=X3)
Physical Surface("Target") = Surface In BoundingBox{
    X3-eps, Y_base-eps, Z1-eps, 
    X3+eps, Y_base+H_tube+eps, Z2+eps
};

// Outlet 2 (X=X4)
Physical Surface("Open") = Surface In BoundingBox{
    X4-eps, Y_base-eps, Z1-eps, 
    X4+eps, Y_base+H_tube+eps, Z2+eps
};

// =============================================================================
// 4. GERAÇÃO E EXPORTAÇÃO
// =============================================================================
Mesh.Algorithm = 8;
Mesh 3;
Save "muffler3D.msh";